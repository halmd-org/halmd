/*
 * Copyright © 2010-2012  Peter Colberg
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

#include <luaponte/class_info.hpp>
#include <stdexcept>

#include <halmd/io/logger.hpp>
#include <halmd/script.hpp>
#include <halmd/utility/filesystem.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/version.h>

using namespace boost;
using namespace std;

HALMD_LUA_API int luaopen_halmd(lua_State* L);

namespace halmd {

script::script()
    // create Lua state, returns raw pointer
  : L(luaL_newstate())
    // wrap Lua state with boost::shared_ptr, and use lua_close as deleter
    // this ensure that any class member variables declared *after* the
    // wrapper L_ will be deconstructed *before* the Lua state,
    // e.g. the Lua script function captured from the HALMD script
  , L_(L, lua_close)
{
    // load Lua standard libraries
    luaL_openlibs(L);
    // load Luabind into Lua interpreter
    load_luaponte();
    // set Lua package path
    package_path();
    // set Lua C package path
    package_cpath();
    // load HALMD Lua C++ wrapper
    luaopen_halmd(L);
    // prepare Lua 5.2 compatible environment
    lua_compat();
}

/**
 * Translate C++ exception into Lua error message
 */
static void translate_exception(lua_State* L, std::exception const& e)
{
    lua_pushstring(L, e.what());
}

/**
 * Load Luabind into Lua interpreter
 */
void script::load_luaponte()
{
    using namespace luaponte;
    // setup global structures and Lua class support
    open(L);
    // print Lua stack trace on error
    set_pcall_callback(&script::traceback);
    // class_info(), class_names()
    bind_class_info(L);
    // translate C++ exception into Lua error message
    register_exception_handler<std::exception>(&translate_exception);
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
 * Set Lua package cpath
 *
 * Append HALMD installation prefix paths to package.cpath.
 */
void script::package_cpath()
{
    // push table "package"
    lua_getglobal(L, "package");
    // push key for rawset
    lua_pushliteral(L, "cpath");
    // push key for rawget
    lua_pushliteral(L, "cpath");
    // get default package.cpath
    lua_rawget(L, -3);

    // absolute path to HALMD build tree
    filesystem::path build_path(HALMD_BINARY_DIR);
    // absolute path to initial current working directory
    filesystem::path initial_path(filesystem::initial_path());

    if (contains_path(build_path, initial_path)) {
        // search for Lua scripts in build tree using relative path
        lua_pushliteral(L, ";" HALMD_BINARY_DIR "/?.so");
    }
    else {
        // search for Lua scripts in installation prefix
        lua_pushliteral(L, ";" HALMD_INSTALL_PREFIX "/lib/?.so");
    }

    // append above literals to default package.cpath
    lua_concat(L, 2);
    // set new package.cpath
    lua_rawset(L, -3);
    // remove table "package"
    lua_pop(L, 1);
}

/**
 * Prepare Lua 5.2 compatible environment
 *
 * http://www.lua.org/manual/5.2/manual.html#8.2
 */
void script::lua_compat()
{
    using namespace luaponte;

#if LUA_VERSION_NUM < 502
    // function unpack was moved into the table library
    globals(L)["table"]["unpack"] = globals(L)["unpack"];
    globals(L)["unpack"] = nil;
#endif

    // function module is deprecated
    globals(L)["module"] = nil;
}

/*
 * Load and execute Lua script
 *
 * If filename is empty, loads from standard input.
 */
void script::dofile(string const& filename)
{
    using namespace luaponte;

    // if filename is NULL, luaL_loadfile loads from standard input
    char const* fn = NULL;
    if (!filename.empty()) {
        fn = filename.c_str();
    }
    // error handler passed to lua_pcall as last argument
    lua_pushcfunction(L, &script::traceback);

    if (luaL_loadfile(L, fn) || lua_pcall(L, 0, 0, 1)) {
        string error(lua_tostring(L, -1));
        lua_pop(L, 1);
        throw runtime_error(error);
    }
}

/**
 * Append traceback to error message on stack
 *
 * @param L Lua state with error message on top of stack
 */
int script::traceback(lua_State* L)
{
    lua_getglobal(L, "debug");
    lua_pushliteral(L, "traceback");
    lua_rawget(L, -2);
    lua_remove(L, -2);
    lua_pushvalue(L, -2);
    lua_remove(L, -3);
    lua_pushnumber(L, 2);
    lua_call(L, 2, 1);
    return 1;
}

} // namespace halmd
