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

#ifndef HALMD_UTILITY_LUA_MAP_CONVERTER_HPP
#define HALMD_UTILITY_LUA_MAP_CONVERTER_HPP

#include <luabind/luabind.hpp>
#include <map>

namespace luabind {

/**
 * Luabind converter for STL map
 */
template <typename Key, typename T>
struct default_converter<std::map<Key, T> >
  : native_converter_base<std::map<Key, T> >
{

    //! compute Lua to C++ conversion score
    static int compute_score(lua_State* L, int index)
    {
        return lua_type(L, index) == LUA_TTABLE ? 0 : -1;
    }

    //! convert from Lua to C++
    std::map<Key, T> from(lua_State* L, int index)
    {
        std::map<Key, T> m;
        for (iterator i(object(from_stack(L, index))), end; i != end; ++i) {
            m[object_cast<Key>(i.key())] = object_cast<T>(*i);
        }
        return m;
    }

    //! convert from C++ to Lua
    void to(lua_State* L, std::map<Key, T> const& m)
    {
        object table = newtable(L);
        typename std::map<Key, T>::const_iterator i, end = m.end();
        for (i = m.begin(); i != end; ++i) {
            // default_converter<T> only invoked with reference wrapper
            table[i->first] = boost::cref(i->second);
        }
        table.push(L);
    }
};

template <typename Key, typename T>
struct default_converter<std::map<Key, T> const&>
  : default_converter<std::map<Key, T> > {};

} // namespace luabind

#endif /* ! HALMD_UTILITY_LUA_MAP_CONVERTER_HPP */
