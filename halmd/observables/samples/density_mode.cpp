/*
 * Copyright © 2011-2013 Felix Höfling
 * Copyright © 2011-2012 Peter Colberg
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

#include <halmd/observables/samples/blocking_scheme.hpp>
#include <halmd/observables/samples/density_mode.hpp>
#include <halmd/utility/lua/lua.hpp>

namespace halmd {
namespace observables {
namespace samples {

void density_mode::luaopen(lua_State* L)
{
    using namespace luaponte;
    module(L)
    [
        namespace_("observables")
        [
            namespace_("samples")
            [
                class_<density_mode>()
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_observables_samples_density_mode(lua_State* L)
{
    density_mode::luaopen(L);
    blocking_scheme<density_mode>::luaopen(L);
    return 0;
}

template class blocking_scheme<density_mode>;

} // namespace samples
} // namespace observables
} // namespace halmd
