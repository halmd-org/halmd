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

#ifndef HALMD_UTILITY_LUA_VECTOR_CONVERTER_HPP
#define HALMD_UTILITY_LUA_VECTOR_CONVERTER_HPP

#include <luabind/luabind.hpp>
#include <vector>

#if LUA_VERSION_NUM < 502
# define luaL_len lua_objlen
#endif

namespace luabind {

/**
 * Luabind converter for STL vector
 */
template <typename T>
struct default_converter<std::vector<T> >
  : native_converter_base<std::vector<T> >
{

    //! compute Lua to C++ conversion score
    static int compute_score(lua_State* L, int index)
    {
        if (lua_type(L, index) != LUA_TTABLE) {
            return -1;
        }
        for (iterator i(object(from_stack(L, index))), end; i != end; ++i) {
            if (!object_cast_nothrow<T>(*i)) {
                return -1;
            }
        }
        return 0;
    }

    //! convert from Lua to C++
    std::vector<T> from(lua_State* L, int index)
    {
        std::size_t len = luaL_len(L, index);
        object table(from_stack(L, index));
        std::vector<T> v;
        v.reserve(len);
        for (std::size_t i = 0; i < len; ++i) {
            v.push_back(object_cast<T>(table[i + 1]));
        }
        return v;
    }

    //! convert from C++ to Lua
    void to(lua_State* L, std::vector<T> const& v)
    {
        object table = newtable(L);
        for (std::size_t i = 0; i < v.size(); ++i) {
            // default_converter<T> only invoked with reference wrapper
            table[i + 1] = boost::cref(v[i]);
        }
        table.push(L);
    }
};

template <typename T>
struct default_converter<std::vector<T> const&>
  : default_converter<std::vector<T> > {};

} // namespace luabind

#if LUA_VERSION_NUM < 502
# undef luaL_len
#endif

#endif /* ! HALMD_UTILITY_LUA_VECTOR_CONVERTER_HPP */
