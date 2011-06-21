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
#include <halmd/utility/filesystem.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/lua/program_options.hpp>
#include <halmd/version.h>

using namespace boost;
using namespace std;

HALMD_LUA_API int luaopen_libhalmd(lua_State* L);

namespace halmd
{

script::script()
  : L(luaL_newstate()) //< create Lua state
{
    luaL_openlibs(L); //< load Lua standard libraries

    package_path(); //< set Lua package path

    load_wrapper(); //< load HALMD Lua C++ wrapper

    load_library(); //< load HALMD Lua library
}

script::~script()
{
    lua_close(L);
}

/**
 * Set Lua package path
 *
 * Append HALMD installation prefix paths to package.path.
 */
void script::package_path()
{
    // push table "package"
    lua_getglobal(L, "package");
    // push key for rawset
    lua_pushliteral(L, "path");
    // push key for rawget
    lua_pushliteral(L, "path");
    // get default package.path
    lua_rawget(L, -3);

    // absolute path to HALMD build tree
    filesystem::path build_path(HALMD_BINARY_DIR);
    // absolute path to initial current working directory
    filesystem::path initial_path(filesystem::initial_path());

    if (contains_path(build_path, initial_path)) {
        // search for Lua scripts in build tree using relative path
        lua_pushliteral(L, ";" HALMD_BINARY_DIR "/lua/?.lua");
        lua_pushliteral(L, ";" HALMD_BINARY_DIR "/lua/?/init.lua");
    }
    else {
        // search for Lua scripts in installation prefix
        lua_pushliteral(L, ";" HALMD_INSTALL_PREFIX "/share/?.lua");
        lua_pushliteral(L, ";" HALMD_INSTALL_PREFIX "/share/?/init.lua");
    }

    // append above literals to default package.path
    lua_concat(L, 3);
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
    using namespace luabind;

    open(L); //< setup global structures and Lua class support

    bind_class_info(L); //< class_info(), class_names()

    luaopen_libhalmd(L);
}

/**
 * Load HALMD Lua library
 */
void script::load_library()
{
    using namespace luabind;

    try {
#ifndef NDEBUG
        scoped_pcall_callback pcall_callback(&traceback);
#endif
        call_function<void>(L, "require", "halmd");
    }
    catch (luabind::error const& e) {
        LOG_ERROR(lua_tostring(L, -1));
        lua_pop(L, 1); //< remove error message
        throw;
    }
}

/*
 * Load and execute Lua script
 */
void script::dofile(string const& file_name)
{
    // error handler passed to lua_pcall as last argument
    lua_pushcfunction(L, &script::traceback);

    if (luaL_loadfile(L, file_name.c_str()) || lua_pcall(L, 0, 0, 1)) {
        LOG_ERROR(lua_tostring(L, -1));
        lua_pop(L, 1); //< remove error message
        throw runtime_error("failed to load Lua script");
    }
}

/**
 * Assemble program options
 */
void script::options(options_parser& parser)
{
    using namespace luabind;

    // retrieve the Lua function before the try-catch block
    // to avoid bogus error message on the Lua stack in case
    // call_function throws an exception
    object options(globals(L)["halmd"]["option"]["get"]);
    try {
#ifndef NDEBUG
        scoped_pcall_callback pcall_callback(&traceback);
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
    using namespace luabind;

    object options(globals(L)["halmd"]["option"]["set"]);
    try {
#ifndef NDEBUG
        scoped_pcall_callback pcall_callback(&traceback);
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
    using namespace luabind;

    // runner is non-template base class to template sampler
    shared_ptr<runner> sampler;
    try {
#ifndef NDEBUG
        scoped_pcall_callback pcall_callback(&traceback);
#endif
        sampler = call_function<shared_ptr<runner> >(L, "halmd");
    }
    catch (luabind::error const& e) {
        LOG_ERROR(lua_tostring(e.state(), -1));
        lua_pop(e.state(), 1); //< remove error message
        throw;
    }

    // Some C++ modules are only needed during the Lua script stage,
    // e.g. the trajectory reader. To make sure these modules are
    // destructed before running the simulation, invoke the Lua
    // garbage collector now.
    lua_gc(L, LUA_GCCOLLECT, 0);

    return sampler;
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
