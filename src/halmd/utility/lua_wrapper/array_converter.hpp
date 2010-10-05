/*
 * Copyright © 2008-2010  Peter Colberg
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

#ifndef HALMD_UTILITY_LUA_WRAPPER_ARRAY_CONVERTER_HPP
#define HALMD_UTILITY_LUA_WRAPPER_ARRAY_CONVERTER_HPP

#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <luabind/luabind.hpp>

namespace luabind
{

/**
 * Luabind converter for Boost array
 */
template <typename T, std::size_t N>
struct default_converter<boost::array<T, N> >
  : native_converter_base<boost::array<T, N> >
{
    //! compute Lua to C++ conversion score
    static int compute_score(lua_State* L, int index)
    {
        return lua_type(L, index) == LUA_TTABLE ? 0 : -1;
    }

    //! convert from Lua to C++
    boost::array<T, N> from(lua_State* L, int index)
    {
        boost::array<T, N> v;
        object table(from_stack(L, index));
        for (std::size_t i = 0; i < v.size(); ++i) {
            v[i] = object_cast<T>(table[i + 1]);
        }
        return v;
    }

    //! convert from C++ to Lua
    void to(lua_State* L, boost::array<T, N> const& array)
    {
        luabind::object table = luabind::newtable(L);
        for (std::size_t i = 0; i < array.size(); ++i) {
            table[i + 1] = array[i];
        }
        table.push(L);
    }
};

template <typename T, std::size_t N>
struct default_converter<boost::array<T, N> const&>
  : default_converter<boost::array<T, N> > {};

/**
 * Luabind converter for 1-dimensional Boost multi_array
 */
template <typename T>
struct default_converter<boost::multi_array<T, 1> >
  : native_converter_base<boost::multi_array<T, 1> >
{
    //! convert from C++ to Lua
    void to(lua_State* L, boost::multi_array<T, 1> const& array)
    {
        luabind::object table = luabind::newtable(L);
        for (std::size_t i = 0; i < array.size(); ++i) {
            table[i + 1] = array[i];
        }
        table.push(L);
    }
};

template <typename T>
struct default_converter<boost::multi_array<T, 1> const&>
  : default_converter<boost::multi_array<T, 1> > {};

} // namespace luabind

#endif /* ! HALMD_UTILITY_LUA_WRAPPER_ARRAY_CONVERTER_HPP */
