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

#include <memory>
#include <stdexcept>

#include <halmd/io/logger.hpp>
#include <halmd/observables/dynamics/correlation.hpp>
#include <halmd/observables/gpu/dynamics/velocity_autocorrelation.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace observables {
namespace gpu {
namespace dynamics {

template <int dimension, typename data_type>
velocity_autocorrelation<dimension, data_type>::velocity_autocorrelation(
    unsigned int blocks
  , unsigned int threads
)
  // allocate block reduction buffers in GPU and page-locked host memory
  : compute_vacf_(blocks, threads)
{
}

template <int dimension, typename data_type>
unsigned int velocity_autocorrelation<dimension, data_type>::defaults::blocks() {
    return 32;
}

template <int dimension, typename data_type>
unsigned int velocity_autocorrelation<dimension, data_type>::defaults::threads() {
    return 256;
}

template <int dimension, typename data_type>
void velocity_autocorrelation<dimension, data_type>::operator() (
    sample_type const& first
  , sample_type const& second
  , accumulator<result_type>& result
)
{
    auto const& velocity1 = first.data();
    auto const& velocity2 = second.data();
    accumulator<result_type> acc = compute_vacf_(
        std::make_tuple(&*velocity1.begin(), &*velocity2.begin())
      , std::make_tuple(&*velocity1.end())
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
void velocity_autocorrelation<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("dynamics")
            [
                class_<velocity_autocorrelation>()

              , def("velocity_autocorrelation", &select_tcf_by_acquire<velocity_autocorrelation>)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_gpu_dynamics_velocity_autocorrelation(lua_State* L)
{
    velocity_autocorrelation<3, float4>::luaopen(L);
    velocity_autocorrelation<2, float4>::luaopen(L);
    observables::dynamics::correlation<velocity_autocorrelation<3, float4> >::luaopen(L);
    observables::dynamics::correlation<velocity_autocorrelation<2, float4> >::luaopen(L);
    return 0;
}

// explicit instantiation
template class velocity_autocorrelation<3, float4>;
template class velocity_autocorrelation<2, float4>;

} // namespace dynamics
} // namespace gpu

namespace dynamics
{

// explicit instantiation
template class correlation<gpu::dynamics::velocity_autocorrelation<3, float4> >;
template class correlation<gpu::dynamics::velocity_autocorrelation<2, float4> >;

} // namespace dynamics
} // namespace observables
} // namespace halmd
