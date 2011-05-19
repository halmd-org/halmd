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

#include <halmd/observables/host/dynamics/velocity_autocorrelation.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace std;

namespace halmd
{
namespace observables { namespace host { namespace dynamics
{

/**
 * Compute velocity autocorrelation of two velocity sample vectors.
 *
 * @param first particles velocities of one species at time t1
 * @param second particles velocities of one species at time t2
 * @returns accumulated velocity autocorrelation
 */
template <int dimension, typename float_type>
typename velocity_autocorrelation<dimension, float_type>::result_type
velocity_autocorrelation<dimension, float_type>::compute(
    sample_vector const& first
  , sample_vector const& second
)
{
    result_type acc;
    typename sample_vector::const_iterator v2, v1, end = first.end();
    for (v2 = first.begin(), v1 = second.begin(); v2 != end; ++v2, ++v1) {
        // accumulate velocity autocorrelation
        acc(inner_prod(*v1, *v2));
    }
    return acc;
}

template <int dimension, typename float_type>
void velocity_autocorrelation<dimension, float_type>::luaopen(lua_State* L)
{
    using namespace luabind;
    static string class_name("velocity_autocorrelation_" + lexical_cast<string>(dimension) + "_");
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("host")
            [
                namespace_("dynamics")
                [
                    class_<velocity_autocorrelation>(class_name.c_str())
                        .def(constructor<>())
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_host_dynamics_velocity_autocorrelation(lua_State* L)
{
#ifndef USE_HOST_SINGLE_PRECISION
    velocity_autocorrelation<3, double>::luaopen(L);
    velocity_autocorrelation<2, double>::luaopen(L);
#else
    velocity_autocorrelation<3, float>::luaopen(L);
    velocity_autocorrelation<2, float>::luaopen(L);
#endif
    return 0;
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class velocity_autocorrelation<3, double>;
template class velocity_autocorrelation<2, double>;
#else
template class velocity_autocorrelation<3, float>;
template class velocity_autocorrelation<2, float>;
#endif

}}} // namespace observables::host::dynamics

} // namespace halmd
