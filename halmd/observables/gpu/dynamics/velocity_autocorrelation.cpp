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

#include <boost/make_shared.hpp>
#include <stdexcept>

#include <halmd/io/logger.hpp>
#include <halmd/observables/dynamics/correlation.hpp>
#include <halmd/observables/gpu/dynamics/velocity_autocorrelation.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace observables {
namespace gpu {
namespace dynamics {

template <int dimension, typename float_type>
velocity_autocorrelation<dimension, float_type>::velocity_autocorrelation(
    unsigned int blocks
  , unsigned int threads
)
  // allocate block reduction buffers in GPU and page-locked host memory
  : compute_vacf_(blocks, threads)
{
}

template <int dimension, typename float_type>
unsigned int velocity_autocorrelation<dimension, float_type>::defaults::blocks() {
    return 32;
}

template <int dimension, typename float_type>
unsigned int velocity_autocorrelation<dimension, float_type>::defaults::threads() {
    return 32 << DEVICE_SCALE;
}

template <int dimension, typename float_type>
typename velocity_autocorrelation<dimension, float_type>::accumulator_type
velocity_autocorrelation<dimension, float_type>::compute(
    sample_type const& first
  , sample_type const& second
)
{
    typename sample_type::velocity_array_type const& velocity1 = first.velocity();
    typename sample_type::velocity_array_type const& velocity2 = second.velocity();
    return compute_vacf_(velocity1.begin(), velocity1.end(), velocity2.begin())();
}

template <typename tcf_type>
static boost::shared_ptr<tcf_type>
select_tcf_by_acquire(boost::function<boost::shared_ptr<typename tcf_type::sample_type const> ()> const&)
{
    return boost::make_shared<tcf_type>();
}

template <int dimension, typename float_type>
void velocity_autocorrelation<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
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
    velocity_autocorrelation<3, float>::luaopen(L);
    velocity_autocorrelation<2, float>::luaopen(L);
    observables::dynamics::correlation<velocity_autocorrelation<3, float> >::luaopen(L);
    observables::dynamics::correlation<velocity_autocorrelation<2, float> >::luaopen(L);
    return 0;
}

// explicit instantiation
template class velocity_autocorrelation<3, float>;
template class velocity_autocorrelation<2, float>;

} // namespace dynamics
} // namespace gpu

namespace dynamics
{

// explicit instantiation
template class correlation<gpu::dynamics::velocity_autocorrelation<3, float> >;
template class correlation<gpu::dynamics::velocity_autocorrelation<2, float> >;

} // namespace dynamics
} // namespace observables
} // namespace halmd
