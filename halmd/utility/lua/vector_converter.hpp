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

namespace luabind
{

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
        return lua_type(L, index) == LUA_TTABLE ? 0 : -1;
    }

    //! convert from Lua to C++
    std::vector<T> from(lua_State* L, int index)
    {
        std::vector<T> v;
        for (iterator i(object(from_stack(L, index))), end; i != end; ++i) {
            v.push_back(object_cast<T>(*i));
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

#endif /* ! HALMD_UTILITY_LUA_VECTOR_CONVERTER_HPP */
