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

#ifndef HALMD_UTILITY_LUA_WRAPPER_ANY_CONVERTER_HPP
#define HALMD_UTILITY_LUA_WRAPPER_ANY_CONVERTER_HPP

#include <boost/any.hpp>
#include <luabind/detail/convert_to_lua.hpp>
#include <luabind/luabind.hpp>
#include <map>
#include <typeinfo>

#include <halmd/utility/demangle.hpp>

// This code is based on the any_converter example of Luabind,
// which demonstrates the conversion of boost::any to Lua.

namespace halmd
{
namespace lua_wrapper { namespace detail
{

typedef void (*any_converter)(lua_State* L, boost::any const&);
typedef std::map<std::type_info const*, any_converter> any_converter_map;

/**
 * singleton registry of boost::any converters
 */
struct any_converters
{
    //! as local static variable to avoid static initialization order fiasco
    static any_converter_map& get()
    {
        static any_converter_map converters;
        return converters;
    }
};

template <typename T>
struct convert_any
{
    //! convert from C++ to Lua
    static void to(lua_State* L, boost::any const& any)
    {
        luabind::detail::convert_to_lua(L, *boost::any_cast<T>(&any));
    }
};

}} // namespace lua_wrapper::detail

/**
 * register a boost::any converter for given type
 */
template <typename T>
void register_any_converter()
{
    using namespace lua_wrapper::detail;
    any_converter_map& converters = any_converters::get();
    converters[&typeid(T)] = convert_any<T>::to;
}

} // namespace halmd

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
    static void to(lua_State* L, boost::any const& any)
    {
        using namespace halmd::lua_wrapper::detail;
        any_converter_map& converters = any_converters::get();
        any_converter_map::const_iterator it = converters.find(&any.type());
        if (it == converters.end()) {
            throw std::runtime_error(
                "unregistered any_converter: " + halmd::demangled_name(any.type())
            );
        }
        it->second(L, any);
    }
};

template <>
struct default_converter<boost::any const&>
  : default_converter<boost::any> {};

} // namespace luabind

#endif /* ! HALMD_UTILITY_LUA_WRAPPER_ANY_CONVERTER_HPP */
