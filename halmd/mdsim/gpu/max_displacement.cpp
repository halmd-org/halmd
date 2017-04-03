/*
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/gpu/max_displacement.hpp>
#include <halmd/utility/lua/lua.hpp>

#include <algorithm>
#include <exception>

namespace halmd {
namespace mdsim {
namespace gpu {

/**
 * construct maximum displacement module
 *
 * @param particle mdsim::gpu::particle instance
 * @param box mdsim::box instance
 */
template <int dimension, typename float_type>
max_displacement<dimension, float_type>::max_displacement(
    std::shared_ptr<particle_type const> particle
  , std::shared_ptr<box_type const> box
)
  // dependency injection
  : particle_(particle)
  , box_(box)
  // select thread-dependent reduction kernel
  , dim_reduce_(64, (64 << DEVICE_SCALE))
  , displacement_impl_(get_displacement_impl(dim_reduce_.threads_per_block()))
  // allocate parameters
  , g_r0_(particle_->nparticle())
  , g_rr_(dim_reduce_.blocks_per_grid())
  , h_rr_(g_rr_.size())
  , displacement_(0)
{
}

template <int dimension, typename float_type>
typename max_displacement_wrapper<dimension>::displacement_impl_type
max_displacement<dimension, float_type>::get_displacement_impl(int threads)
{
    switch (threads) {
      case 512:
        return max_displacement_wrapper<dimension>::kernel.displacement_impl[0];
      case 256:
        return max_displacement_wrapper<dimension>::kernel.displacement_impl[1];
      case 128:
        return max_displacement_wrapper<dimension>::kernel.displacement_impl[2];
      case 64:
        return max_displacement_wrapper<dimension>::kernel.displacement_impl[3];
      case 32:
        return max_displacement_wrapper<dimension>::kernel.displacement_impl[4];
      default:
        throw std::logic_error("invalid reduction thread count");
    }
}

/**
 * Zero maximum squared displacement
 */
template <int dimension, typename float_type>
void max_displacement<dimension, float_type>::zero()
{
    cache<position_array_type> const& position_cache = particle_->position();
    cuda::vector<float4> const& position = read_cache(position_cache);

    LOG_TRACE("zero maximum squared displacement");

    scoped_timer_type timer(runtime_.zero);
    cuda::copy(position.begin(), position.begin() + particle_->nparticle(), g_r0_.begin());
    displacement_ = 0;
    position_cache_ = position_cache;
}

/**
 * Compute maximum squared displacement
 */
template <int dimension, typename float_type>
float_type max_displacement<dimension, float_type>::compute()
{
    cache<position_array_type> const& position_cache = particle_->position();

    if (position_cache != position_cache_) {
        position_array_type const& position = read_cache(position_cache);

        LOG_TRACE("compute maximum squared displacement");

        scoped_timer_type timer(runtime_.compute);
        try {
            cuda::configure(
                dim_reduce_.grid
              , dim_reduce_.block
              , dim_reduce_.threads_per_block() * sizeof(float)
            );
            displacement_impl_(
                position.data()
              , g_r0_
              , g_rr_
              , particle_->nparticle()
              , static_cast<vector_type>(box_->length())
            );
            cuda::copy(g_rr_, h_rr_);
        }
        catch (cuda::error const&) {
            LOG_ERROR("failed to reduce squared particle displacements on GPU");
            throw;
        }
        displacement_ = std::sqrt(*std::max_element(h_rr_.begin(), h_rr_.end()));
        position_cache_ = position_cache;
    }
    return displacement_;
}

template <int dimension, typename float_type>
void max_displacement<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            class_<max_displacement>()
                .scope
                [
                    class_<runtime>("runtime")
                        .def_readonly("zero", &runtime::zero)
                        .def_readonly("compute", &runtime::compute)
                ]
                .def_readonly("runtime", &max_displacement::runtime_)
          , def("max_displacement", &std::make_shared<max_displacement
                  , std::shared_ptr<particle_type const>
                  , std::shared_ptr<box_type const>
              >)
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_max_displacement(lua_State* L)
{
#ifdef USE_GPU_SINGLE_PRECISION
    max_displacement<3, float>::luaopen(L);
    max_displacement<2, float>::luaopen(L);
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
    max_displacement<3, dsfloat>::luaopen(L);
    max_displacement<2, dsfloat>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifdef USE_GPU_SINGLE_PRECISION
template class max_displacement<3, float>;
template class max_displacement<2, float>;
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
template class max_displacement<3, dsfloat>;
template class max_displacement<2, dsfloat>;
#endif

} // namespace mdsim
} // namespace gpu
} // namespace halmd
