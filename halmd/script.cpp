/*
 * Copyright © 2008-2010  Peter Colberg and Felix Höfling
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

#include <boost/bind.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/script.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>
#include <halmd/utility/lua_wrapper/ublas.hpp>
#include <halmd/utility/lua_wrapper/variables_map.hpp>
#include <halmd/version.h>

using namespace boost;
using namespace std;

namespace halmd
{

script::script()
  : L_(luaL_newstate(), lua_close) //< create Lua state
{
    lua_State* L = get_pointer(L_); //< get raw pointer for Lua C API

    using namespace luabind;

    luaL_openlibs(L); //< load Lua standard libraries

    set_pcall_callback(&traceback); //< set pcall error handler

    package_path(); //< set Lua package path

    load_wrapper(); //< load HALMD Lua C++ wrapper

    load_library(); //< load HALMD Lua library
}

/**
 * Set Lua package path
 */
void script::package_path()
{
    lua_State* L = get_pointer(L_); //< get raw pointer for Lua C API

    using namespace luabind;

    string path;
    path.append( HALMD_BINARY_DIR "/lua/?.lua" ";" );
    path.append( HALMD_BINARY_DIR "/lua/?/init.lua" ";" );
    path.append( HALMD_SOURCE_DIR "/lua/?.lua" ";" );
    path.append( HALMD_SOURCE_DIR "/lua/?/init.lua" ";" );
    path.append( HALMD_INSTALL_PREFIX "/share/?.lua" ";" );
    path.append( HALMD_INSTALL_PREFIX "/share/?/init.lua" ";" );
    path.append( HALMD_INSTALL_PREFIX "/lib/?.lua" ";" );
    path.append( HALMD_INSTALL_PREFIX "/lib/?/init.lua" ";" );
    path.append( object_cast<string>(globals(L)["package"]["path"]) );
    globals(L)["package"]["path"] = path;
}

/**
 * Load HALMD Lua wrapper
 *
 * Register C++ classes with Lua.
 */
void script::load_wrapper()
{
    lua_State* L = get_pointer(L_); //< get raw pointer for Lua C API

    using namespace luabind;

    open(L); //< setup global structures and Lua class support

    lua_wrapper::open(L); //< register HALMD Lua wrappers

    lua_wrapper::ublas::luaopen(L);
}

/**
 * Load HALMD Lua library
 */
void script::load_library()
{
    lua_State* L = get_pointer(L_); //< get raw pointer for Lua C API

    using namespace luabind;

    try {
        call_function<void>(L, "require", "halmd");
    }
    catch (luabind::error const& e) {
        LOG_ERROR(lua_tostring(L, -1));
        lua_pop(L, 1); //< remove error message
        throw;
    }
}

/**
 * Assemble program options
 */
po::options_description script::options()
{
    lua_State* L = get_pointer(L_); //< get raw pointer for Lua C API

    using namespace luabind;

    po::options_description desc;

    // retrieve the Lua function before the try-catch block
    // to avoid bogus error message on the Lua stack in case
    // call_function throws an exception
    object options(globals(L)["halmd"]["modules"]["options"]);
    try {
        call_function<void>(options, ref(desc));
    }
    catch (luabind::error const& e) {
        LOG_ERROR(lua_tostring(e.state(), -1));
        lua_pop(e.state(), 1); //< remove error message
        throw;
    }

    return desc; //< avoid variable in main() at the expense of return-by-value
}

/**
 * Set parsed command line options
 */
void script::parsed(po::variables_map const& vm)
{
    lua_State* L = get_pointer(L_); //< get raw pointer for Lua C API

    using namespace luabind;

    object options(globals(L)["halmd"]["modules"]["parsed"]);
    try {
        call_function<void>(options, cref(vm));
    }
    catch (luabind::error const& e) {
        LOG_ERROR(lua_tostring(e.state(), -1));
        lua_pop(e.state(), 1); //< remove error message
        throw;
    }
}

/**
 * Run simulation
 */
void script::run()
{
    lua_State* L = get_pointer(L_); //< get raw pointer for Lua C API

    using namespace luabind;

    try {
        call_function<void>(L, "run");
    }
    catch (luabind::error const& e) {
        LOG_ERROR(lua_tostring(e.state(), -1));
        lua_pop(e.state(), 1); //< remove error message
        throw;
    }
}

/**
 * Append traceback to error message on stack
 *
 * @param L Lua state with error message on top of stack
 */
int script::traceback(lua_State* L)
{
    lua_pushliteral(L, "\n");
    lua_getfield(L, LUA_GLOBALSINDEX, "debug");
    lua_pushliteral(L, "traceback");
    lua_rawget(L, -2);
    lua_remove(L, -2);
    lua_call(L, 0, 1);
    lua_concat(L, 3);
    return 1;
}

} // namespace halmd
