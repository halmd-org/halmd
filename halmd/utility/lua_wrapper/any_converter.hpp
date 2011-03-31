/*
 * Copyright © 2010-2011  Peter Colberg
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

#ifndef HALMD_UTILITY_LUA_WRAPPER_ANY_CONVERTER_HPP
#define HALMD_UTILITY_LUA_WRAPPER_ANY_CONVERTER_HPP

#include <boost/any.hpp>
#include <luabind/luabind.hpp>

namespace luabind
{

/**
 * Luabind converter for boost::any
 */
template <>
struct default_converter<boost::any>
  : native_converter_base<boost::any>
{
    //! convert from C++ to Lua
    static void to(lua_State* L, boost::any const& any);
};

template <>
struct default_converter<boost::any const&>
  : default_converter<boost::any> {};

} // namespace luabind

#endif /* ! HALMD_UTILITY_LUA_WRAPPER_ANY_CONVERTER_HPP */
