/*
 * Copyright © 2008-2010  Peter Colberg
 * Copyright © 2013       Nicolas Höft
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

#include <algorithm>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <cmath>

#include <halmd/algorithm/gpu/radix_sort.hpp>
#include <halmd/mdsim/gpu/sorts/hilbert.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace halmd::algorithm::gpu;

namespace halmd {
namespace mdsim {
namespace gpu {
namespace sorts {

template <int dimension, typename float_type>
hilbert<dimension, float_type>::hilbert(
    std::shared_ptr<particle_type> particle
  , std::shared_ptr<box_type const> box
  , std::shared_ptr<logger> logger
)
  // dependency injection
  : particle_(particle)
  , box_(box)
  , logger_(logger)
{
    // FIXME set Hilbert space-filling curve recursion depth
    float max_length = *std::max_element(box_->length().begin(), box_->length().end());
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
    LOG_TRACE("order particles after Hilbert space-filling curve");
    {
        scoped_timer_type timer(runtime_.order);
        cuda::vector<unsigned int> g_index(particle_->nparticle());
        g_index.reserve(particle_->dim().threads());
        {

            cuda::vector<unsigned int> g_map(particle_->nparticle());
            g_map.reserve(particle_->dim().threads());
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
    position_array_type const& position = read_cache(particle_->position());

    scoped_timer_type timer(runtime_.map);

    int blockSize = wrapper_type::kernel.map.max_block_size();
    if (!blockSize) blockSize = particle_->dim().block.x;
    cuda::configure(particle_->array_size() / blockSize, blockSize);
    wrapper_type::kernel.map(
        position.data()
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
    int blockSize = wrapper_type::kernel.gen_index.max_block_size();
    if (!blockSize) blockSize = particle_->dim().block.x;
    cuda::configure(particle_->array_size() / blockSize, blockSize);
    wrapper_type::kernel.gen_index(g_index);
    radix_sort(g_map.begin(), g_map.end(), g_index.begin());
}

template <typename sort_type>
static std::function<void ()>
wrap_order(std::shared_ptr<sort_type> self)
{
    return [=]() {
        self->order();
    };
}

template <int dimension, typename float_type>
void hilbert<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("mdsim")
        [
            namespace_("sorts")
            [
                class_<hilbert>()
                    .property("order", &wrap_order<hilbert>)
                    .def("on_order", &hilbert::on_order)
                    .scope
                    [
                        class_<runtime>("runtime")
                            .def_readonly("order", &runtime::order)
                            .def_readonly("map", &runtime::map)
                    ]
                    .def_readonly("runtime", &hilbert::runtime_)
              , def("hilbert", &std::make_shared<hilbert
                    , std::shared_ptr<particle_type>
                    , std::shared_ptr<box_type const>
                    , std::shared_ptr<logger>
                >)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_mdsim_gpu_sorts_hilbert(lua_State* L)
{
#ifdef USE_GPU_SINGLE_PRECISION
    hilbert<3, float>::luaopen(L);
    hilbert<2, float>::luaopen(L);
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
    hilbert<3, dsfloat>::luaopen(L);
    hilbert<2, dsfloat>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifdef USE_GPU_SINGLE_PRECISION
template class hilbert<3, float>;
template class hilbert<2, float>;
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
template class hilbert<3, dsfloat>;
template class hilbert<2, dsfloat>;
#endif

} // namespace sorts
} // namespace gpu
} // namespace mdsim
} // namespace halmd
