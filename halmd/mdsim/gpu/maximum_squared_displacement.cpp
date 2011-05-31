/*
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <algorithm>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/range/iterator_range.hpp>
#include <exception>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/gpu/maximum_squared_displacement.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace mdsim { namespace gpu
{

/**
 * construct maximum displacement module
 *
 * @param particle mdsim::gpu::particle instance
 * @param box mdsim::box instance
 */
template <int dimension, typename float_type>
maximum_squared_displacement<dimension, float_type>::maximum_squared_displacement(
    shared_ptr<particle_type const> particle
  , shared_ptr<box_type const> box
)
  // dependency injection
  : particle_(particle)
  , box_(box)
  // select thread-dependent reduction kernel
  , dim_reduce_(64, (64 << DEVICE_SCALE))
  , displacement_impl_(get_displacement_impl(dim_reduce_.threads_per_block()))
  // allocate parameters
  , g_r0_(particle_->nbox)
  , g_rr_(dim_reduce_.blocks_per_grid())
  , h_rr_(g_rr_.size())
{
    try {
        cuda::copy(particle_->nbox, get_maximum_squared_displacement_kernel<dimension>().nbox);
        cuda::copy(static_cast<vector_type>(box_->length()), get_maximum_squared_displacement_kernel<dimension>().box_length);
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to copy parameters to device symbols");
        throw;
    }
}

/**
 * register module runtime accumulators
 */
template <int dimension, typename float_type>
void maximum_squared_displacement<dimension, float_type>::register_runtimes(profiler_type& profiler)
{
    profiler.register_runtime(runtime_.zero, "zero", "zero maximum displacement");
    profiler.register_runtime(runtime_.compute, "compute", "compute maximum displacement");
}

template <int dimension, typename float_type>
typename maximum_squared_displacement_wrapper<dimension>::displacement_impl_type
maximum_squared_displacement<dimension, float_type>::get_displacement_impl(int threads)
{
    switch (threads) {
      case 512:
        return maximum_squared_displacement_wrapper<dimension>::kernel.displacement_impl[0];
      case 256:
        return maximum_squared_displacement_wrapper<dimension>::kernel.displacement_impl[1];
      case 128:
        return maximum_squared_displacement_wrapper<dimension>::kernel.displacement_impl[2];
      case 64:
        return maximum_squared_displacement_wrapper<dimension>::kernel.displacement_impl[3];
      case 32:
        return maximum_squared_displacement_wrapper<dimension>::kernel.displacement_impl[4];
      default:
        throw std::logic_error("invalid reduction thread count");
    }
}

/**
 * Zero maximum squared displacement
 */
template <int dimension, typename float_type>
void maximum_squared_displacement<dimension, float_type>::zero()
{
    LOG_TRACE("zero maximum squared displacement");

    scoped_timer<timer> timer_(runtime_.zero);
    cuda::copy(particle_->g_r, g_r0_);
}

/**
 * Compute maximum squared displacement
 */
template <int dimension, typename float_type>
float_type maximum_squared_displacement<dimension, float_type>::compute()
{
    LOG_TRACE("compute maximum squared displacement");

    scoped_timer<timer> timer_(runtime_.compute);
    try {
        cuda::configure(dim_reduce_.grid, dim_reduce_.block, dim_reduce_.threads_per_block() * sizeof(float));
        displacement_impl_(particle_->g_r, g_r0_, g_rr_);
        cuda::copy(g_rr_, h_rr_);
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to reduce squared particle displacements on GPU");
        throw;
    }
    return *max_element(h_rr_.begin(), h_rr_.end());
}

template <int dimension, typename float_type>
void maximum_squared_displacement<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("maximum_squared_displacement_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("gpu")
            [
                class_<maximum_squared_displacement, shared_ptr<maximum_squared_displacement> >(class_name.c_str())
                    .def(constructor<
                        shared_ptr<particle_type const>
                      , shared_ptr<box_type const>
                    >())
                    .def("register_runtimes", &maximum_squared_displacement::register_runtimes)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_maximum_squared_displacement(lua_State* L)
{
    maximum_squared_displacement<3, float>::luaopen(L);
    maximum_squared_displacement<2, float>::luaopen(L);
    return 0;
}

// explicit instantiation
template class maximum_squared_displacement<3, float>;
template class maximum_squared_displacement<2, float>;

}} // namespace mdsim::gpu

} // namespace halmd
