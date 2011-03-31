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

#include <luabind/luabind.hpp>
#include <luabind/detail/convert_to_lua.hpp>
#include <typeinfo>

#include <halmd/config.hpp>
#include <halmd/utility/demangle.hpp>
#include <halmd/utility/lua/any_converter.hpp>
#include <halmd/utility/lua/array_converter.hpp>
#include <halmd/utility/lua/long_long_converter.hpp>
#include <halmd/utility/lua/program_options.hpp>
#include <halmd/utility/lua/vector_converter.hpp>

using namespace boost;
using namespace std;

// This code is based on the any_converter example of Luabind,
// which demonstrates the conversion of boost::any to Lua.

namespace luabind
{

//! convert from C++ to Lua
void default_converter<any>::to(lua_State* L, any const& value)
{
    typedef void (*converter)(lua_State* L, any const&);

    // get table of any converter functions from Lua registry,
    // using address of default_converter<any>::to as unique key
    lua_pushlightuserdata(L, (void*)&to);
    lua_gettable(L, LUA_REGISTRYINDEX);

    // lookup any converter function by typeid address
    lua_pushlightuserdata(L, (void*)&value.type());
    lua_gettable(L, -2);
    converter convert = (converter)lua_touserdata(L, -1);
    lua_pop(L, 2);
    if (!convert) {
        throw runtime_error("unregistered any converter: " + halmd::demangled_name(value.type()));
    }

    // convert any to value type and push onto stack
    convert(L, value);
}

} // namespace luabind

namespace halmd
{

template <typename T>
struct any_converter
{
    //! convert from C++ to Lua
    static void to(lua_State* L, any const& any)
    {
        luabind::detail::convert_to_lua(L, *any_cast<T>(&any));
    }
};

template <>
struct any_converter<void>
{
    //! convert from C++ to Lua
    static void to(lua_State* L, any const&)
    {
        lua_pushnil(L);
    }
};

/**
 * register an any converter for given type
 *
 * expects table on top of stack
 */
template <typename T>
static void register_any_converter(lua_State* L)
{
    lua_pushlightuserdata(L, (void*)&typeid(T));
    lua_pushlightuserdata(L, (void*)&any_converter<T>::to);
    lua_settable(L, -3);
}

HALMD_LUA_API int luaopen_libhalmd_any_converter(lua_State* L)
{
    // create table of any converter functions in Lua registry,
    // using address of default_converter<any>::to as unique key
    lua_newtable(L);
    lua_pushlightuserdata(L, (void*)&luabind::default_converter<any>::to);
    lua_pushvalue(L, -2);
    lua_settable(L, LUA_REGISTRYINDEX);

    register_any_converter<void>(L); //< empty any

    register_any_converter<bool>(L);
    register_any_converter<char>(L);
    register_any_converter<signed char>(L);
    register_any_converter<unsigned char>(L);
    register_any_converter<signed short>(L);
    register_any_converter<unsigned short>(L);
    register_any_converter<signed int>(L);
    register_any_converter<unsigned int>(L);
    register_any_converter<signed long>(L);
    register_any_converter<unsigned long>(L);
    register_any_converter<signed long long>(L);
    register_any_converter<unsigned long long>(L);
    register_any_converter<float>(L);
    register_any_converter<double>(L);
    register_any_converter<long double>(L);
    register_any_converter<string>(L);
    register_any_converter<char const*>(L);

    register_any_converter<vector<bool> >(L);
    register_any_converter<vector<char> >(L);
    register_any_converter<vector<signed char> >(L);
    register_any_converter<vector<unsigned char> >(L);
    register_any_converter<vector<signed short> >(L);
    register_any_converter<vector<unsigned short> >(L);
    register_any_converter<vector<signed int> >(L);
    register_any_converter<vector<unsigned int> >(L);
    register_any_converter<vector<signed long> >(L);
    register_any_converter<vector<unsigned long> >(L);
    register_any_converter<vector<signed long long> >(L);
    register_any_converter<vector<unsigned long long> >(L);
    register_any_converter<vector<float> >(L);
    register_any_converter<vector<double> >(L);
    register_any_converter<vector<long double> >(L);
    register_any_converter<vector<string> >(L);
    register_any_converter<vector<char const*> >(L);

    register_any_converter<multi_array<bool, 1> >(L);
    register_any_converter<multi_array<char, 1> >(L);
    register_any_converter<multi_array<signed char, 1> >(L);
    register_any_converter<multi_array<unsigned char, 1> >(L);
    register_any_converter<multi_array<signed short, 1> >(L);
    register_any_converter<multi_array<unsigned short, 1> >(L);
    register_any_converter<multi_array<signed int, 1> >(L);
    register_any_converter<multi_array<unsigned int, 1> >(L);
    register_any_converter<multi_array<signed long, 1> >(L);
    register_any_converter<multi_array<unsigned long, 1> >(L);
    register_any_converter<multi_array<signed long long, 1> >(L);
    register_any_converter<multi_array<unsigned long long, 1> >(L);
    register_any_converter<multi_array<float, 1> >(L);
    register_any_converter<multi_array<double, 1> >(L);
    register_any_converter<multi_array<long double, 1> >(L);
    register_any_converter<multi_array<string, 1> >(L);
    register_any_converter<multi_array<char const*, 1> >(L);

    register_any_converter<program_options::variables_map>(L);

    lua_pop(L, 1);
    return 0;
}

} // namespace halmd
