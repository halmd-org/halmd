/*
 * Copyright © 2013 Felix Höfling
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

#include <halmd/io/logger.hpp>
#include <halmd/observables/dynamics/correlation.hpp>
#include <halmd/observables/dynamics/intermediate_scattering_function.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace std;

namespace halmd {
namespace observables {
namespace dynamics {

template <int dimension>
void intermediate_scattering_function<dimension>::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("dynamics")
            [
                class_<intermediate_scattering_function>()

              , def("intermediate_scattering_function", &make_shared<intermediate_scattering_function
                  , shared_ptr<wavevector_type const>
                  , double
                >)
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_dynamics_intermediate_scattering_function(lua_State* L)
{
    intermediate_scattering_function<3>::luaopen(L);
    intermediate_scattering_function<2>::luaopen(L);
    observables::dynamics::correlation<intermediate_scattering_function<3>>::luaopen(L);
    observables::dynamics::correlation<intermediate_scattering_function<2>>::luaopen(L);
    return 0;
}

// explicit instantiation
template class intermediate_scattering_function<3>;
template class intermediate_scattering_function<2>;

// explicit instantiation
template class correlation<intermediate_scattering_function<3>>;
template class correlation<intermediate_scattering_function<2>>;

} // namespace dynamics
} // namespace observables
} // namespace halmd
