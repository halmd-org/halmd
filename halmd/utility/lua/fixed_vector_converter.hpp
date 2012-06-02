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

#ifndef HALMD_UTILITY_LUA_FIXED_VECTOR_CONVERTER_HPP
#define HALMD_UTILITY_LUA_FIXED_VECTOR_CONVERTER_HPP

#include <algorithm>
#include <boost/bind.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <luaponte/luaponte.hpp>

#include <halmd/config.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>

#if LUA_VERSION_NUM < 502
# define luaL_len lua_objlen
#endif

namespace luaponte {

/**
 * Luabind converter for fixed-size algebraic vector
 */
template <typename T, std::size_t N>
struct default_converter<halmd::fixed_vector<T, N> >
  : native_converter_base<halmd::fixed_vector<T, N> >
{

    //! compute Lua to C++ conversion score
    static int compute_score(lua_State* L, int index)
    {
        return lua_type(L, index) == LUA_TTABLE && luaL_len(L, index) == N ? 0 : -1;
    }

    //! convert from Lua to C++
    halmd::fixed_vector<T, N> from(lua_State* L, int index)
    {
        halmd::fixed_vector<T, N> v;
        object table(from_stack(L, index));
        for (std::size_t i = 0; i < v.size(); ++i) {
            v[i] = object_cast<T>(table[i + 1]);
        }
        return v;
    }

    //! convert from C++ to Lua
    void to(lua_State* L, halmd::fixed_vector<T, N> const& v)
    {
        object table = newtable(L);
        for (std::size_t i = 0; i < v.size(); ++i) {
            // default_converter<T> only invoked with reference wrapper
            table[i + 1] = boost::cref(v[i]);
        }
        table.push(L);
    }
};

template <typename T, std::size_t N>
struct default_converter<halmd::fixed_vector<T, N> const&>
  : default_converter<halmd::fixed_vector<T, N> > {};

#ifndef HALMD_NO_CXX11
template <typename T, std::size_t N>
struct default_converter<halmd::fixed_vector<T, N>&&>
  : default_converter<halmd::fixed_vector<T, N> > {};
#endif

} // namespace luaponte

#if LUA_VERSION_NUM < 502
# undef luaL_len
#endif

#endif /* ! HALMD_UTILITY_LUA_FIXED_VECTOR_CONVERTER_HPP */
