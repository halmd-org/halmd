/*
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

#ifndef HALMD_UTILITY_LUA_CACHE_CONVERTER_HPP
#define HALMD_UTILITY_LUA_CACHE_CONVERTER_HPP

#include <boost/ref.hpp>
#include <luaponte/luaponte.hpp>

#if LUA_VERSION_NUM < 502
# define luaL_len lua_objlen
#endif

namespace halmd {

// forward declaration
template <typename T>
class cache;
template <typename T>
class cache_proxy;

} // namespace halmd

namespace luaponte {

/**
 * Lua value converter for halmd::cache
 */
template <typename T>
struct default_converter<halmd::cache<T> >
  : native_converter_base<halmd::cache<T> >
{
    //! convert from C++ to Lua
    void to(lua_State* L, halmd::cache<T> const& cache)
    {
        halmd::cache_proxy<T const> value = cache;
        object(L, boost::cref(*value)).push(L);
    }
};

template <typename T>
struct default_converter<halmd::cache<T> const&>
  : default_converter<halmd::cache<T> > {};

template <typename T>
struct default_converter<halmd::cache<T>&&>
  : default_converter<halmd::cache<T> > {};

template <typename T>
struct default_converter<halmd::cache<T>&>
  : default_converter<halmd::cache<T> > {};

} // namespace luaponte

#if LUA_VERSION_NUM < 502
# undef luaL_len
#endif

#endif /* ! HALMD_UTILITY_LUA_CACHE_CONVERTER_HPP */
