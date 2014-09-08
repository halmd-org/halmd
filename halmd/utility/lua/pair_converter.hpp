/*
 * Copyright © 2014 Felix Höfling
 * Copyright © 2012 Peter Colberg
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

#ifndef HALMD_UTILITY_LUA_PAIR_CONVERTER_HPP
#define HALMD_UTILITY_LUA_PAIR_CONVERTER_HPP

#include <luaponte/luaponte.hpp>

#include <utility>

#if LUA_VERSION_NUM < 502
# define luaL_len lua_objlen
#endif

namespace luaponte {

/**
 * Lua table converter for std::pair
 */
template <typename T, typename U>
struct default_converter<std::pair<T, U> >
  : native_converter_base<std::pair<T, U> >
{
    //! compute Lua to C++ conversion score
    static int compute_score(lua_State* L, int index)
    {
        return luaL_len(L, index) == 2 ? 0 : -1;
    }

    //! convert from Lua to C++
    std::pair<T, U> from(lua_State* L, int index)
    {
        object table(from_stack(L, index));
        return std::make_pair(object_cast<T>(table[1]), object_cast<U>(table[2]));
    }

    //! convert from C++ to Lua
    void to(lua_State* L, std::pair<T, U> const& v)
    {
        object table = newtable(L);
        table[1] = boost::cref(v.first);
        table[2] = boost::cref(v.second);
        table.push(L);
    }
};

template <typename T, typename U>
struct default_converter<std::pair<T, U> const&>
  : default_converter<std::pair<T, U> > {};

template <typename T, typename U>
struct default_converter<std::pair<T, U>&&>
  : default_converter<std::pair<T, U> > {};

template <typename T, typename U>
struct default_converter<std::pair<T, U>&>
  : default_converter<std::pair<T, U> > {};

} // namespace luaponte

#if LUA_VERSION_NUM < 502
# undef luaL_len
#endif

#endif /* ! HALMD_UTILITY_LUA_PAIR_CONVERTER_HPP */
