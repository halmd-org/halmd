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

namespace halmd
{

template <typename Signature>
static void luaopen_signal_proxy(lua_State* L, char const* class_name)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        class_<signal_proxy<Signature> >(class_name)
            .scope
            [
                class_<typename signal_proxy<Signature>::connection>("connection")

              , class_<typename signal_proxy<Signature>::slot_function_type>("slot_function_type")
                    .def("__call", &signal_proxy<Signature>::slot_function_type::operator())
            ]
            .def("connect", &signal_proxy<Signature>::connect)
            .def("disconnect", &signal_proxy<Signature>::disconnect)
    ];
}

HALMD_LUA_API int luaopen_libhalmd_utility_lua_signal(lua_State* L)
{
    luaopen_signal_proxy<void ()>(L, "signal_proxy__void__");
    luaopen_signal_proxy<void (uint64_t)>(L, "signal_proxy__uint64_t__");
    return 0;
}

} // namespace halmd
