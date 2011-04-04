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

#ifndef HALMD_UTILITY_LUA_OPTIONAL_CONVERTER_HPP
#define HALMD_UTILITY_LUA_OPTIONAL_CONVERTER_HPP

#include <boost/optional.hpp>
#include <luabind/luabind.hpp>

namespace luabind
{

/**
 * Luabind converter for Boost.Optional
 */
template <typename T>
struct default_converter<boost::optional<T> >
  : native_converter_base<boost::optional<T> >
{
    //! compute Lua to C++ conversion score
    static int compute_score(lua_State* L, int index)
    {
        return luabind::object_cast_nothrow<T>(from_stack(L, index)) ? 0 : -1;
    }

    //! convert from Lua to C++
    boost::optional<T> from(lua_State* L, int index)
    {
        return luabind::object_cast<T>(from_stack(L, index));
    }

    //! convert from C++ to Lua
    void to(lua_State* L, boost::optional<T> const& value)
    {
        if (value) {
            luabind::detail::convert_to_lua(L, *value);
        }
        else {
            lua_pushnil(L);
        }
    }
};

template <typename T>
struct default_converter<boost::optional<T> const&>
  : default_converter<boost::optional<T> > {};

} // namespace luabind

#endif /* ! HALMD_UTILITY_LUA_OPTIONAL_CONVERTER_HPP */
