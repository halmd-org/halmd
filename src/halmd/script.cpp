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
#include <halmd/version.h>

using namespace boost;
using namespace std;

namespace halmd
{

script::script()
  : L_(luaL_newstate(), lua_close) //< create Lua state
{
    lua_State* L = get_pointer(L_); //< get raw pointer for Lua C API

    luaL_openlibs(L); //< load Lua standard libraries

    load_wrapper(); //< load HALMD Lua C++ wrapper

    load_library(); //< load HALMD Lua library
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
}

/**
 * Load HALMD Lua library
 */
void script::load_library()
{
    lua_State* L = get_pointer(L_); //< get raw pointer for Lua C API

    using namespace luabind;

    string path;
    path.append( HALMD_BINARY_DIR "/lib/?.lua" ";" );
    path.append( HALMD_BINARY_DIR "/lib/?/init.lua" ";" );
    path.append( HALMD_SOURCE_DIR "/lib/?.lua" ";" );
    path.append( HALMD_SOURCE_DIR "/lib/?/init.lua" ";" );
    path.append( HALMD_INSTALL_PREFIX "/share/?.lua" ";" );
    path.append( HALMD_INSTALL_PREFIX "/share/?/init.lua" ";" );
    path.append( HALMD_INSTALL_PREFIX "/lib/?.lua" ";" );
    path.append( HALMD_INSTALL_PREFIX "/lib/?/init.lua" ";" );
    path.append( object_cast<string>(globals(L)["package"]["path"]) );
    globals(L)["package"]["path"] = path;

    try {
        call_function<void>(L, "require", "halmd");
    }
    catch (luabind::error const& e) {
        LOG_ERROR("[Lua] " << lua_tostring(L, -1));
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

    po::options_description desc("Program Options");

    // retrieve the Lua function before the try-catch block
    // to avoid bogus error message on the Lua stack in case
    // call_function throws an exception
    object options(globals(L)["halmd"]["modules"]["options"]);
    try {
        call_function<void>(options, ref(desc));
    }
    catch (luabind::error const& e) {
        LOG_ERROR("[Lua] " << lua_tostring(e.state(), -1));
        lua_pop(e.state(), 1); //< remove error message
        throw;
    }

    return desc; //< avoid variable in main() at the expense of return-by-value
}

/**
 * Set parsed command line options
 */
void script::options(po::variables_map const& vm)
{
    lua_State* L = get_pointer(L_); //< get raw pointer for Lua C API

    using namespace luabind;

    object options(globals(L)["halmd"]["options"]);
    try {
        call_function<void>(options, cref(vm));
    }
    catch (luabind::error const& e) {
        LOG_ERROR("[Lua] " << lua_tostring(e.state(), -1));
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

    object run(globals(L)["halmd"]["run"]);
    try {
        call_function<void>(run);
    }
    catch (luabind::error const& e) {
        LOG_ERROR("[Lua] " << lua_tostring(e.state(), -1));
        lua_pop(e.state(), 1); //< remove error message
        throw;
    }
}

} // namespace halmd
