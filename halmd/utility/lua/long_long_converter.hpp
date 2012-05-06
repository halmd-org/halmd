/*
 * Copyright © 2010  Peter Colberg
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

#ifndef HALMD_UTILITY_LUA_LONG_LONG_CONVERTER_HPP
#define HALMD_UTILITY_LUA_LONG_LONG_CONVERTER_HPP

#include <luabind/luabind.hpp>

#include <halmd/config.hpp>

namespace luabind {

/**
 * Luabind converter for long long integer
 */
template <>
struct default_converter<long long>
  : native_converter_base<long long>
{
    //! compute Lua to C++ conversion score
    static int compute_score(lua_State* L, int index)
    {
        return lua_type(L, index) == LUA_TNUMBER ? 0 : -1;
    }

    //! convert from Lua to C++
    long long from(lua_State* L, int index)
    {
        return static_cast<long long>(lua_tonumber(L, index));
    }

    //! convert from C++ to Lua
    void to(lua_State* L, long long value)
    {
        lua_pushnumber(L, static_cast<lua_Number>(value));
    }
};

template <>
struct default_converter<long long const&>
  : default_converter<long long> {};

#ifndef HALMD_NO_CXX11
template <>
struct default_converter<long long&&>
  : default_converter<long long> {};
#endif

/**
 * Luabind converter for unsigned long long integer
 */
template <>
struct default_converter<unsigned long long>
  : native_converter_base<unsigned long long>
{
    //! compute Lua to C++ conversion score
    static int compute_score(lua_State* L, int index)
    {
        return lua_type(L, index) == LUA_TNUMBER ? 0 : -1;
    }

    //! convert from Lua to C++
    unsigned long long from(lua_State* L, int index)
    {
        return static_cast<unsigned long long>(lua_tonumber(L, index));
    }

    //! convert from C++ to Lua
    void to(lua_State* L, unsigned long long value)
    {
        lua_pushnumber(L, static_cast<lua_Number>(value));
    }
};

template <>
struct default_converter<unsigned long long const&>
  : default_converter<unsigned long long> {};

#ifndef HALMD_NO_CXX11
template <>
struct default_converter<unsigned long long&&>
  : default_converter<unsigned long long> {};
#endif

} // namespace luabind

#endif /* ! HALMD_UTILITY_LUA_LONG_LONG_CONVERTER_HPP */
