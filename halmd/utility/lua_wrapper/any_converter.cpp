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

#include <luabind/luabind.hpp>
#include <halmd/utility/lua_wrapper/any_converter.hpp>
#include <halmd/utility/lua_wrapper/array_converter.hpp>
#include <halmd/utility/lua_wrapper/long_long_converter.hpp>
#include <halmd/utility/lua_wrapper/vector_converter.hpp>

namespace halmd
{

static __attribute__((constructor)) void register_any_converters()
{
    register_any_converter<void>(); //< empty boost::any

    register_any_converter<bool>();
    register_any_converter<char>();
    register_any_converter<signed char>();
    register_any_converter<unsigned char>();
    register_any_converter<signed short>();
    register_any_converter<unsigned short>();
    register_any_converter<signed int>();
    register_any_converter<unsigned int>();
    register_any_converter<signed long>();
    register_any_converter<unsigned long>();
    register_any_converter<signed long long>();
    register_any_converter<unsigned long long>();
    register_any_converter<float>();
    register_any_converter<double>();
    register_any_converter<long double>();
    register_any_converter<std::string>();
    register_any_converter<char const*>();

    register_any_converter<std::vector<bool> >();
    register_any_converter<std::vector<char> >();
    register_any_converter<std::vector<signed char> >();
    register_any_converter<std::vector<unsigned char> >();
    register_any_converter<std::vector<signed short> >();
    register_any_converter<std::vector<unsigned short> >();
    register_any_converter<std::vector<signed int> >();
    register_any_converter<std::vector<unsigned int> >();
    register_any_converter<std::vector<signed long> >();
    register_any_converter<std::vector<unsigned long> >();
    register_any_converter<std::vector<signed long long> >();
    register_any_converter<std::vector<unsigned long long> >();
    register_any_converter<std::vector<float> >();
    register_any_converter<std::vector<double> >();
    register_any_converter<std::vector<long double> >();
    register_any_converter<std::vector<std::string> >();
    register_any_converter<std::vector<char const*> >();

    register_any_converter<boost::multi_array<bool, 1> >();
    register_any_converter<boost::multi_array<char, 1> >();
    register_any_converter<boost::multi_array<signed char, 1> >();
    register_any_converter<boost::multi_array<unsigned char, 1> >();
    register_any_converter<boost::multi_array<signed short, 1> >();
    register_any_converter<boost::multi_array<unsigned short, 1> >();
    register_any_converter<boost::multi_array<signed int, 1> >();
    register_any_converter<boost::multi_array<unsigned int, 1> >();
    register_any_converter<boost::multi_array<signed long, 1> >();
    register_any_converter<boost::multi_array<unsigned long, 1> >();
    register_any_converter<boost::multi_array<signed long long, 1> >();
    register_any_converter<boost::multi_array<unsigned long long, 1> >();
    register_any_converter<boost::multi_array<float, 1> >();
    register_any_converter<boost::multi_array<double, 1> >();
    register_any_converter<boost::multi_array<long double, 1> >();
    register_any_converter<boost::multi_array<std::string, 1> >();
    register_any_converter<boost::multi_array<char const*, 1> >();
}

} // namespace halmd
