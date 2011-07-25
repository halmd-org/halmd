/*
 * Copyright © 2011  Peter Colberg
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

#include <luabind/luabind.hpp>
#include <stdint.h>

#include <halmd/config.hpp>
#include <halmd/utility/signal.hpp>

namespace halmd {

HALMD_LUA_API int luaopen_libhalmd_utility_lua_signal(lua_State* L)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        class_<signal<void ()> >("signal__void__")
            .scope
            [
                class_<signal<void ()>::slot_function_type>("slot_function_type")
                    .def("__call", &signal<void ()>::slot_function_type::operator())
            ]
    ];
    return 0;
}

} // namespace halmd
