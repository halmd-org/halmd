/*
 * Copyright © 2008-2010  Peter Colberg
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
#include <cmath>

#include <halmd/algorithm/gpu/radix_sort.hpp>
#include <halmd/mdsim/gpu/sorts/hilbert.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace halmd::algorithm::gpu;
using namespace std;

namespace halmd {
namespace mdsim {
namespace gpu {
namespace sorts {

template <int dimension, typename float_type>
hilbert<dimension, float_type>::hilbert(
    boost::shared_ptr<particle_type> particle
  , boost::shared_ptr<box_type const> box
  , boost::shared_ptr<logger_type> logger
)
  // dependency injection
  : particle_(particle)
  , box_(box)
  , logger_(logger)
{
    // FIXME set Hilbert space-filling curve recursion depth
    float_type max_length = *max_element(box_->length().begin(), box_->length().end());
    depth_ = static_cast<unsigned int>(std::ceil(std::log(max_length) / M_LN2));
    // 32-bit integer for 2D/3D Hilbert code allows a maximum of 16/10 levels
    depth_ = std::min((dimension == 3) ? 10U : 16U, depth_);

    LOG("vertex recursion depth: " << depth_);

    try {
        cuda::copy(depth_, wrapper_type::kernel.depth);
    }
    catch (cuda::error const&) {
        LOG_ERROR("failed to copy parameters to device");
        throw;
    }
}

/**
 * Order particles after Hilbert space-filling curve
 */
template <int dimension, typename float_type>
void hilbert<dimension, float_type>::order()
{
    LOG_TRACE("order particles");
    {
        cuda::vector<unsigned int> g_index(particle_->nparticle());
        g_index.reserve(particle_->position().capacity());
        scoped_timer_type timer(runtime_.order);
        {

            cuda::vector<unsigned int> g_map(particle_->nparticle());
            g_map.reserve(particle_->position().capacity());
            this->map(g_map);
            this->permutation(g_map, g_index);
        }
        particle_->rearrange(g_index);
    }
    on_order_();
}

/**
 * map particles to Hilbert curve
 */
template <int dimension, typename float_type>
void hilbert<dimension, float_type>::map(cuda::vector<unsigned int>& g_map)
{
    scoped_timer_type timer(runtime_.map);
    cuda::configure(particle_->dim.grid, particle_->dim.block);
    wrapper_type::kernel.map(
        particle_->position()
      , g_map
      , static_cast<vector_type>(box_->length())
    );
}

/**
 * generate permutation
 */
template <int dimension, typename float_type>
void hilbert<dimension, float_type>::permutation(cuda::vector<unsigned int>& g_map, cuda::vector<unsigned int>& g_index)
{
    cuda::configure(particle_->dim.grid, particle_->dim.block);
    wrapper_type::kernel.gen_index(g_index);
    radix_sort<unsigned int> sort(particle_->nparticle(), particle_->dim.threads_per_block());
    sort(g_map, g_index);
}

template <int dimension, typename float_type>
static char const* module_name_wrapper(hilbert<dimension, float_type> const&)
{
    return hilbert<dimension, float_type>::module_name();
}

template <typename sort_type>
static boost::function<void ()>
wrap_order(boost::shared_ptr<sort_type> sort)
{
    return boost::bind(&sort_type::order, sort);
}

template <int dimension, typename float_type>
void hilbert<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("hilbert_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("gpu")
            [
                namespace_("sorts")
                [
                    class_<hilbert, boost::shared_ptr<hilbert> >(class_name.c_str())
                        .def(constructor<
                            boost::shared_ptr<particle_type>
                          , boost::shared_ptr<box_type const>
                          , boost::shared_ptr<logger_type>
                        >())
                        .property("module_name", &module_name_wrapper<dimension, float_type>)
                        .property("order", &wrap_order<hilbert>)
                        .def("on_order", &hilbert::on_order)
                        .scope
                        [
                            class_<runtime>("runtime")
                                .def_readonly("order", &runtime::order)
                                .def_readonly("map", &runtime::map)
                        ]
                        .def_readonly("runtime", &hilbert::runtime_)
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_sorts_hilbert(lua_State* L)
{
    hilbert<3, float>::luaopen(L);
    hilbert<2, float>::luaopen(L);
    return 0;
}

// explicit instantiation
template class hilbert<3, float>;
template class hilbert<2, float>;

} // namespace mdsim
} // namespace gpu
} // namespace sorts
} // namespace halmd
