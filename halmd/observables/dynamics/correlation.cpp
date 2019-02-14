/*
 * Copyright © 2011 Felix Höfling
 * Copyright © 2011-2012 Peter Colberg
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

#include <halmd/observables/dynamics/correlation.hpp>

namespace halmd {
namespace observables {
namespace dynamics {

void correlation_base::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L, "libhalmd")
    [
        namespace_("observables")
        [
            namespace_("dynamics")
            [
                class_<correlation_base>()
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_dynamics_correlation(lua_State* L)
{
    correlation_base::luaopen(L);
    return 0;
}

} // namespace dynamics
} // namespace observables
} // namespace halmd
