/*
 * Copyright © 2014 Felix Höfling
 * Copyright © 2012 Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#ifndef HALMD_UTILITY_LUA_RAW_ARRAY_CONVERTER_HPP
#define HALMD_UTILITY_LUA_RAW_ARRAY_CONVERTER_HPP

#include <luaponte/luaponte.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/utility/demangle.hpp>

#if LUA_VERSION_NUM < 502
# define luaL_len lua_objlen
#endif

namespace halmd {

// forward declaration
template <typename T>
class raw_array;

} // namespace halmd

namespace luaponte {

/**
 * Luabind converter for STL vector
 */
template <typename T>
struct default_converter<halmd::raw_array<T> >
  : native_converter_base<halmd::raw_array<T> >
{

    /**
     * Compute Lua to C++ conversion score.
     */
    static int compute_score(lua_State* L, int index)
    {
        return lua_type(L, index) == LUA_TTABLE ? 0 : -1;
    }

    /**
     * Convert from Lua to C++.
     */
    halmd::raw_array<T> from(lua_State* L, int index)
    {
        std::size_t len = luaL_len(L, index);
        object table(from_stack(L, index));
        LOG_TRACE("convert Lua table of size " << len << " to halmd::raw_array<" << demangled_name<T>() << ">");
        halmd::raw_array<T> v(len);
        for (std::size_t i = 0; i < len; ++i) {
            v[i] = object_cast<T>(table[i + 1]);
        }
        return v;
    }

    /**
     * Convert from C++ to Lua.
     */
    void to(lua_State* L, halmd::raw_array<T> const& v)
    {
        LOG_TRACE("convert halmd::raw_array<" << demangled_name<T>() << "> of size " << v.size() << " to Lua table");
        object table = newtable(L);
        std::size_t i = 1;
        for (auto const& x : v) {
            // default_converter<T> only invoked with reference wrapper
            table[i++] = boost::cref(x);
        }
        table.push(L);
    }
};

template <typename T>
struct default_converter<halmd::raw_array<T> const&>
  : default_converter<halmd::raw_array<T> > {};

template <typename T>
struct default_converter<halmd::raw_array<T>&&>
  : default_converter<halmd::raw_array<T> > {};

template <typename T>
struct default_converter<halmd::raw_array<T>&>
  : default_converter<halmd::raw_array<T> > {};

} // namespace luaponte

#if LUA_VERSION_NUM < 502
# undef luaL_len
#endif

#endif /* ! HALMD_UTILITY_LUA_RAW_ARRAY_CONVERTER_HPP */
