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

#include <lua.hpp>

#include <halmd/config.hpp>

HALMD_LUA_API int luaopen_libhalmd_any_converter(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_h5(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_po(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_signal(lua_State* L);
HALMD_LUA_API int luaopen_libhalmd_ublas(lua_State* L);

HALMD_LUA_API int luaopen_libhalmd(lua_State* L)
{
    luaopen_libhalmd_any_converter(L);
    luaopen_libhalmd_h5(L);
    luaopen_libhalmd_po(L);
    luaopen_libhalmd_signal(L);
    luaopen_libhalmd_ublas(L);
    return 0;
}
