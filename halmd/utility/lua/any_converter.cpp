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
#include <luabind/typeid.hpp>

#include <halmd/config.hpp>
#include <halmd/utility/demangle.hpp>
#include <halmd/utility/lua/any_converter.hpp>
#include <halmd/utility/lua/array_converter.hpp>
#include <halmd/utility/lua/long_long_converter.hpp>
#include <halmd/utility/lua/program_options.hpp>
#include <halmd/utility/lua/vector_converter.hpp>

// This code is based on the any_converter example of Luabind,
// which demonstrates the conversion of boost::any to Lua.
//
// The implementation below, however, does not use &typeid(T)
// to store a type T as a key in the any converter map. This is
// not portable (e.g. to AIX), as there is no guarantee that
// &typeid(T) == &typeid(T) for identical type T.
// Instead we use the Luabind class type_id, which uses the
// well-defined method type_info::before to compare types.
//
// The any converter map is stored inside the Lua interpreter,
// not as a global variable, to allow unit testing with multiple
// Lua interpreter instantiations in a single test executable.

using namespace boost;
using namespace std;

namespace halmd {

class any_converter
  : public map<luabind::type_id, void (*)(lua_State* L, any const&)>
{
private:
    template <typename T>
    struct convert
    {
        //! convert from C++ to Lua
        static void to(lua_State* L, any const& any)
        {
            luabind::detail::convert_to_lua(L, *any_cast<T>(&any));
        }
    };

public:
    template <typename T>
    void register_type()
    {
        this->insert(make_pair(luabind::type_id(typeid(T)), &convert<T>::to));
    }
};

template <>
struct any_converter::convert<void>
{
    //! convert from C++ to Lua
    static void to(lua_State* L, any const&)
    {
        lua_pushnil(L);
    }
};

HALMD_LUA_API int luaopen_libhalmd_utility_lua_any_converter(lua_State* L)
{
    using namespace luabind;

    module(L, "libhalmd")
    [
        class_<any_converter>("any_converter")
    ];

    registry(L)["libhalmd_any_converter"] = any_converter();

    any_converter* conv = object_cast<any_converter*>(registry(L)["libhalmd_any_converter"]);

    conv->register_type<void>(); //< empty any

    conv->register_type<bool>();
    conv->register_type<char>();
    conv->register_type<signed char>();
    conv->register_type<unsigned char>();
    conv->register_type<signed short>();
    conv->register_type<unsigned short>();
    conv->register_type<signed int>();
    conv->register_type<unsigned int>();
    conv->register_type<signed long>();
    conv->register_type<unsigned long>();
    conv->register_type<signed long long>();
    conv->register_type<unsigned long long>();
    conv->register_type<float>();
    conv->register_type<double>();
    conv->register_type<long double>();
    conv->register_type<string>();
    conv->register_type<char const*>();

    conv->register_type<vector<bool> >();
    conv->register_type<vector<char> >();
    conv->register_type<vector<signed char> >();
    conv->register_type<vector<unsigned char> >();
    conv->register_type<vector<signed short> >();
    conv->register_type<vector<unsigned short> >();
    conv->register_type<vector<signed int> >();
    conv->register_type<vector<unsigned int> >();
    conv->register_type<vector<signed long> >();
    conv->register_type<vector<unsigned long> >();
    conv->register_type<vector<signed long long> >();
    conv->register_type<vector<unsigned long long> >();
    conv->register_type<vector<float> >();
    conv->register_type<vector<double> >();
    conv->register_type<vector<long double> >();
    conv->register_type<vector<string> >();
    conv->register_type<vector<char const*> >();

    conv->register_type<multi_array<bool, 1> >();
    conv->register_type<multi_array<char, 1> >();
    conv->register_type<multi_array<signed char, 1> >();
    conv->register_type<multi_array<unsigned char, 1> >();
    conv->register_type<multi_array<signed short, 1> >();
    conv->register_type<multi_array<unsigned short, 1> >();
    conv->register_type<multi_array<signed int, 1> >();
    conv->register_type<multi_array<unsigned int, 1> >();
    conv->register_type<multi_array<signed long, 1> >();
    conv->register_type<multi_array<unsigned long, 1> >();
    conv->register_type<multi_array<signed long long, 1> >();
    conv->register_type<multi_array<unsigned long long, 1> >();
    conv->register_type<multi_array<float, 1> >();
    conv->register_type<multi_array<double, 1> >();
    conv->register_type<multi_array<long double, 1> >();
    conv->register_type<multi_array<string, 1> >();
    conv->register_type<multi_array<char const*, 1> >();

    conv->register_type<program_options::variables_map>();

    return 0;
}

} // namespace halmd

namespace luabind {

//! convert from C++ to Lua
void default_converter<any>::to(lua_State* L, any const& value)
{
    using namespace halmd;
    any_converter const* conv = object_cast<any_converter const*>(registry(L)["libhalmd_any_converter"]);
    if (!conv) {
        throw runtime_error("no registered any converters");
    }
    any_converter::const_iterator it = conv->find(value.type());
    if (it == conv->end()) {
        throw runtime_error("unregistered any converter: " + demangled_name(value.type()));
    }
    it->second(L, value);
}

} // namespace luabind
