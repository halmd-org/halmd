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

#include <memory>
#include <stdexcept>

#include <halmd/io/logger.hpp>
#include <halmd/observables/dynamics/correlation.hpp>
#include <halmd/observables/gpu/dynamics/mean_square_displacement.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace observables {
namespace gpu {
namespace dynamics {

template <int dimension, typename float_type>
mean_square_displacement<dimension, float_type>::mean_square_displacement(
    unsigned int blocks
  , unsigned int threads
)
  // allocate block reduction buffers in GPU and page-locked host memory
  : compute_msd_(blocks, threads)
{
}

template <int dimension, typename float_type>
unsigned int mean_square_displacement<dimension, float_type>::defaults::blocks() {
    return 32;
}

template <int dimension, typename float_type>
unsigned int mean_square_displacement<dimension, float_type>::defaults::threads() {
    return 32 << DEVICE_SCALE;
}

template <int dimension, typename float_type>
void mean_square_displacement<dimension, float_type>::operator() (
    sample_type const& first
  , sample_type const& second
  , accumulator<result_type>& result
)
{
    typename sample_type::position_array_type const& position1 = first.position();
    typename sample_type::position_array_type const& position2 = second.position();
    accumulator<result_type> acc = compute_msd_(
        std::make_tuple(&*position1.begin(), &*position2.begin())
      , std::make_tuple(&*position1.end())
    )();
    result(acc);
}

template <typename tcf_type>
static std::shared_ptr<tcf_type>
select_tcf_by_acquire(std::function<std::shared_ptr<typename tcf_type::sample_type const> ()> const&)
{
    return std::make_shared<tcf_type>();
}

template <int dimension, typename float_type>
void mean_square_displacement<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("dynamics")
            [
                class_<mean_square_displacement>()

              , def("mean_square_displacement", &select_tcf_by_acquire<mean_square_displacement>)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_gpu_dynamics_mean_square_displacement(lua_State* L)
{
    mean_square_displacement<3, float>::luaopen(L);
    mean_square_displacement<2, float>::luaopen(L);
    observables::dynamics::correlation<mean_square_displacement<3, float> >::luaopen(L);
    observables::dynamics::correlation<mean_square_displacement<2, float> >::luaopen(L);
    return 0;
}

// explicit instantiation
template class mean_square_displacement<3, float>;
template class mean_square_displacement<2, float>;

} // namespace dynamics
} // namespace gpu

namespace dynamics
{

// explicit instantiation
template class correlation<gpu::dynamics::mean_square_displacement<3, float> >;
template class correlation<gpu::dynamics::mean_square_displacement<2, float> >;

} // namespace dynamics
} // namespace observables
} // namespace halmd
