/*
 * Copyright © 2010  Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#ifndef HALMD_TEST_TOOLS_LUA_HPP
#define HALMD_TEST_TOOLS_LUA_HPP

#include <lua.hpp>
#include <luaponte/luaponte.hpp>
#include <memory>
#include <ostream>

#define LUA_WARN( str )         BOOST_WARN_MESSAGE( lua_test_fixture::dostring( str ), lua_test_fixture::error(lua_test_fixture::L) )
#define LUA_CHECK( str )        BOOST_CHECK_MESSAGE( lua_test_fixture::dostring( str ), lua_test_fixture::error(lua_test_fixture::L) )
#define LUA_REQUIRE( str )      BOOST_REQUIRE_MESSAGE( lua_test_fixture::dostring( str ), lua_test_fixture::error(lua_test_fixture::L) )

/**
 * Lua test fixture.
 *
 * Setup Lua state with standard libraries and Luabind.
 */
class lua_test_fixture
{
private:
    std::shared_ptr<lua_State> L_;

public:
    lua_State* const L;

    lua_test_fixture();
    bool dostring(std::string const& str);

    struct error
    {
    public:
        lua_State* const L;

        error(lua_State* L) : L(L) {}
    };
};

extern std::ostream& operator<<(std::ostream& os, lua_test_fixture::error const&);

#endif /* ! HALMD_TEST_TOOLS_LUA_HPP */
