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

#define BOOST_TEST_MODULE map_converter
#include <boost/test/unit_test.hpp>

#include "test/utility/lua_wrapper/lua_setup.hpp"

using namespace boost;
using namespace std;

/**
 * Test conversion of std::map from and to Lua.
 */

std::map<std::string, int> static_map_to_lua()
{
    std::map<std::string, int> m;
    m["foo"] = 5u;
    m["bar"] = 42;
    return m;
}

void lua_to_static_map(std::map<std::string, int> const& m)
{
    BOOST_CHECK( m.find("foo")->second == 5 );
    BOOST_CHECK( m.find("bar")->second == 42 );
    BOOST_CHECK( m.find("foobar")->second == 43 );
}

/**
 * Test std::map converter for value type with build-in Lua converter
 */
BOOST_FIXTURE_TEST_CASE( static_map, lua_setup )
{
    using namespace luabind;

    module(L)
    [
        def("static_map_to_lua", &static_map_to_lua)
      , def("lua_to_static_map", &lua_to_static_map)
    ];

    LUA_REQUIRE( "m = assert(static_map_to_lua())" );
    LUA_CHECK( "assert(m.foo == 5)" );
    LUA_CHECK( "assert(m.bar == 42)" );
    LUA_CHECK( "assert(m.foobar == nil)" );
    LUA_CHECK( "m.foobar = 43" );
    LUA_CHECK( "lua_to_static_map(m)" );
}

any any_to_lua() { return 42; }

std::map<std::string, boost::any> any_map_to_lua()
{
    std::map<std::string, boost::any> m;
    m["foo"] = 5u;
    m["bar"] = 42;
    return m;
}

/**
 * Test std::map converter for value type with custom Lua converter
 */
BOOST_FIXTURE_TEST_CASE( any_map, lua_setup )
{
    using namespace luabind;

    module(L)
    [
        def("any_to_lua", &any_to_lua)
      , def("any_map_to_lua", &any_map_to_lua)
      , def("lua_to_static_map", &lua_to_static_map)
    ];

    LUA_REQUIRE( "assert(any_to_lua() == 42)" );
    LUA_REQUIRE( "m = assert(any_map_to_lua())" );
    LUA_CHECK( "assert(m.foo == 5)" );
    LUA_CHECK( "assert(m.bar == 42)" );
    LUA_CHECK( "assert(m.foobar == nil)" );
    LUA_CHECK( "m.foobar = 43" );
    LUA_CHECK( "lua_to_static_map(m)" );
}

static __attribute__((constructor)) void register_any_converters()
{
    using namespace halmd;
    register_any_converter<unsigned int>();
    register_any_converter<int>();
}
