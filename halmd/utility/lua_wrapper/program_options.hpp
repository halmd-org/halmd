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

#ifndef HALMD_UTILITY_LUA_WRAPPER_PROGRAM_OPTIONS_HPP
#define HALMD_UTILITY_LUA_WRAPPER_PROGRAM_OPTIONS_HPP

#include <boost/program_options/variables_map.hpp>
#include <luabind/luabind.hpp>

#include <halmd/utility/lua_wrapper/any_converter.hpp>

namespace halmd
{
namespace lua_wrapper { namespace program_options
{

int luaopen(lua_State* L);

}} // namespace lua_wrapper::program_options

} // namespace halmd

namespace luabind
{

/**
 * Luabind converter for Boost.Program_options variables_map
 */
template <>
struct default_converter<boost::program_options::variables_map>
  : native_converter_base<boost::program_options::variables_map>
{
    //! convert from C++ to Lua
    void to(lua_State* L, boost::program_options::variables_map const& vm)
    {
        luabind::object table = luabind::newtable(L);
        boost::program_options::variables_map::const_iterator it, end = vm.end();
        for (it = vm.begin(); it != end; ++it) {
            // convert boost::any to Lua value
            table[it->first] = boost::cref(it->second.value());
        }
        table.push(L);
    }
};

template <>
struct default_converter<boost::program_options::variables_map const&>
  : default_converter<boost::program_options::variables_map> {};

} // namespace luabind

#endif /* ! HALMD_UTILITY_LUA_WRAPPER_PROGRAM_OPTIONS_HPP */
