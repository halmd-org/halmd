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

#include "test/utility/lua_wrapper/lua_setup.hpp"

/**
 * Initialize Lua state.
 *
 * Load Lua standard libraries and Luabind global structures.
 */
lua_setup::lua_setup()
  : L_(luaL_newstate(), lua_close) // call lua_close upon destruction
  , L(get_pointer(L_))
{
    luaL_openlibs(L);

    luabind::open(L);
}

static int pcall_handler(lua_State* L)
{
    return 1;
}

/**
 * This is a modified version of dostring in the Luabind unit tests
 * adapted to Boost Test. Instead of throwing the Lua error message as
 * an exception, it returns true in case of success and false in case
 * of failure. This value is passed as a predicate to BOOST_*_CHECK.
 */
bool lua_setup::dostring(std::string const& str)
{
    lua_pushcclosure(L, &pcall_handler, 0);

    if (luaL_loadbuffer(L, str.c_str(), str.length(), str.c_str()))
    {
        return false;
    }
    if (lua_pcall(L, 0, 0, -2))
    {
        return false;
    }

    lua_pop(L, 1);

    return true;
}

/**
 * This function retrieves an error message from the Lua stack.
 * It is called in case of failure with BOOST_*_CHECK_MESSAGE.
 */
std::ostream& operator<<(std::ostream& os, lua_setup::error const& e)
{
    os << lua_tostring(e.L, -1);

    lua_pop(e.L, 2);

    return os;
}
