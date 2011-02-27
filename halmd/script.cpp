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
#include <luabind/class_info.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/script.hpp>
#include <halmd/utility/lua_wrapper/error.hpp>
#include <halmd/utility/lua_wrapper/hdf5.hpp>
#include <halmd/utility/lua_wrapper/lua_wrapper.hpp>
#include <halmd/utility/lua_wrapper/program_options.hpp>
#include <halmd/utility/lua_wrapper/ublas.hpp>
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

    package_path(); //< set Lua package path

    load_wrapper(); //< load HALMD Lua C++ wrapper

    load_library(); //< load HALMD Lua library
}

/**
 * Set Lua package path
 *
 * Append HALMD installation prefix paths to package.path.
 */
void script::package_path()
{
    lua_State* L = get_pointer(L_); //< get raw pointer for Lua C API

    // push table "package"
    lua_getglobal(L, "package");
    // push key for rawset
    lua_pushliteral(L, "path");
    // push key for rawget
    lua_pushliteral(L, "path");
    // get default package.path
    lua_rawget(L, -3);

    // search for Lua scripts in build tree using relative path
    lua_pushliteral(L, ";" "../lua/?.lua");
    lua_pushliteral(L, ";" "../lua/?/init.lua");

    // search for Lua scripts in installation prefix
    lua_pushliteral(L, ";" HALMD_INSTALL_PREFIX "/share/?.lua");
    lua_pushliteral(L, ";" HALMD_INSTALL_PREFIX "/share/?/init.lua");
    lua_pushliteral(L, ";" HALMD_INSTALL_PREFIX "/lib/?.lua");
    lua_pushliteral(L, ";" HALMD_INSTALL_PREFIX "/lib/?/init.lua");

    // append above literals to default package.path
    lua_concat(L, 7);
    // set new package.path
    lua_rawset(L, -3);
    // remove table "package"
    lua_pop(L, 1);
}

/**
 * Load HALMD Lua wrapper
 *
 * Register C++ classes with Lua.
 */
void script::load_wrapper()
{
    lua_State* L = get_pointer(L_); //< get raw pointer for Lua C API

    luabind::open(L); //< setup global structures and Lua class support

    luabind::bind_class_info(L); //< class_info(), class_names()

    lua_wrapper::open(L); //< register HALMD Lua wrappers

    lua_wrapper::hdf5::luaopen(L);
    lua_wrapper::program_options::luaopen(L);
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
#ifndef NDEBUG
        lua_wrapper::scoped_pcall_callback pcall_callback(&traceback);
#endif
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
void script::options(options_parser& parser)
{
    lua_State* L = get_pointer(L_); //< get raw pointer for Lua C API

    using namespace luabind;

    // retrieve the Lua function before the try-catch block
    // to avoid bogus error message on the Lua stack in case
    // call_function throws an exception
    object options(globals(L)["halmd"]["modules"]["options"]);
    try {
#ifndef NDEBUG
        lua_wrapper::scoped_pcall_callback pcall_callback(&traceback);
#endif
        call_function<void>(options, ref(parser));
    }
    catch (luabind::error const& e) {
        LOG_ERROR(lua_tostring(e.state(), -1));
        lua_pop(e.state(), 1); //< remove error message
        throw;
    }
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
#ifndef NDEBUG
        lua_wrapper::scoped_pcall_callback pcall_callback(&traceback);
#endif
        call_function<void>(options, vm);
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
shared_ptr<runner> script::run()
{
    lua_State* L = get_pointer(L_); //< get raw pointer for Lua C API

    using namespace luabind;

    try {
#ifndef NDEBUG
        lua_wrapper::scoped_pcall_callback pcall_callback(&traceback);
#endif
        call_function<void>(L, "run");
    }
    catch (luabind::error const& e) {
        LOG_ERROR(lua_tostring(e.state(), -1));
        lua_pop(e.state(), 1); //< remove error message
        throw;
    }

    object sampler(globals(L)["halmd"]["sampler"]);

    // downcast from template class sampler to base class runner
    return call_function<shared_ptr<runner> >(sampler);
}

/**
 * Append traceback to error message on stack
 *
 * @param L Lua state with error message on top of stack
 */
int script::traceback(lua_State* L)
{
    lua_pushliteral(L, "\n");
    lua_getglobal(L, "debug");
    lua_pushliteral(L, "traceback");
    lua_rawget(L, -2);
    lua_remove(L, -2);
    lua_call(L, 0, 1);
    lua_concat(L, 3);
    return 1;
}

} // namespace halmd
