/*
 * Copyright © 2015      Felix Höfling
 * Copyright © 2010-2012 Peter Colberg
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

namespace fs = boost::filesystem;

HALMD_LUA_API int luaopen_halmd(lua_State* L);

namespace halmd {

script::script()
    // create Lua state, returns raw pointer
  : L(luaL_newstate())
    // wrap Lua state with std::shared_ptr, and use lua_close as deleter
    // this ensure that any class member variables declared *after* the
    // wrapper L_ will be deconstructed *before* the Lua state,
    // e.g. the Lua script function captured from the HALMD script
  , L_(L, lua_close)
{
    // load Lua standard libraries
    luaL_openlibs(L);
    // load Luabind into Lua interpreter
    load_luaponte();

    // set Lua package path and C package path
    //
    // We use the HALMD build tree if the initial current working directory is inside
    // the build tree, otherwise the installation tree is used.
    if (contains_path(fs::path(HALMD_BINARY_DIR), fs::path(fs::initial_path()))) {
        // search for Lua scripts in build tree using relative path
        prepend_package_path(HALMD_BINARY_DIR "/lua");
        prepend_package_cpath(HALMD_BINARY_DIR);
    }
    else {
        // search for Lua scripts in installation prefix
        prepend_package_path(HALMD_INSTALL_PREFIX "/share/halmd/lua");
        prepend_package_cpath(HALMD_INSTALL_PREFIX "/lib");
    }

    // load HALMD Lua C++ wrapper
    luaopen_halmd(L);
}

/**
 * Translate C++ exception into Lua error message
 */
static void translate_exception(lua_State* L, std::exception const& e)
{
    lua_pushstring(L, e.what());
}

/**
 * Translate Luaponte exception into Lua error message.
 *
 * Luaponte pushes the error message onto the stack.
 */
static void translate_lua_error(lua_State* L, luaponte::error const&)
{
    lua_pushliteral(L, "[Lua] ");
    lua_pushvalue(L, -2);
    lua_concat(L, 2);
    lua_remove(L, -2);
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
    register_exception_handler<luaponte::error>(&translate_lua_error);
}

/**
 * prepend Lua's package.path by 'path'
 */
void script::prepend_package_path(std::string const& path)
{
    // push table "package"
    lua_getglobal(L, "package");
    // push key for rawset
    lua_pushliteral(L, "path");

    // add Lua scripts in given path and subdirectories
    std::string s = path + "/?.lua;";
    lua_pushlstring(L, s.c_str(), s.size());

    s = path + "/?/init.lua;";
    lua_pushlstring(L, s.c_str(), s.size());

    // push key for rawget
    lua_pushliteral(L, "path");
    // get default package.path
    lua_rawget(L, -5);

    // append the above 2 strings to default package.path
    lua_concat(L, 3);
    // set new package.path
    lua_rawset(L, -3);
    // remove table "package"
    lua_pop(L, 1);
}

/**
 * prepend Lua's package.cpath by 'path'
 */
void script::prepend_package_cpath(std::string const& path)
{
    // push table "package"
    lua_getglobal(L, "package");
    // push key for rawset
    lua_pushliteral(L, "cpath");

    // add shared libraries in given path
    std::string s = path + "/?.so;";
    lua_pushlstring(L, s.c_str(), s.size());

    // push key for rawget
    lua_pushliteral(L, "cpath");
    // get default package.path
    lua_rawget(L, -4);

    // append the above string to default package.path
    lua_concat(L, 2);
    // set new package.path
    lua_rawset(L, -3);
    // remove table "package"
    lua_pop(L, 1);
}

/*
 * Load and execute Lua script
 *
 * If filename is empty, loads from standard input.
 */
void script::dofile(std::string const& filename)
{
    using namespace luaponte;

    // if filename is NULL, luaL_loadfile loads from standard input
    char const* fn = NULL;
    if (!filename.empty()) {
        fn = filename.c_str();
        // search for Lua scripts relative to script directory
        prepend_package_path(fs::path(filename).parent_path().string());
    }

    // error handler passed to lua_pcall as last argument
    lua_pushcfunction(L, &script::traceback);

    if (luaL_loadfile(L, fn) || lua_pcall(L, 0, 0, 1)) {
        std::string error(lua_tostring(L, -1));
        lua_pop(L, 1);
        throw std::runtime_error(error);
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
