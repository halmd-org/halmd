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
#include <halmd/utility/filesystem.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/lua/program_options.hpp>
#include <halmd/version.h>

using namespace boost;
using namespace std;

HALMD_LUA_API int luaopen_libhalmd(lua_State* L);

namespace halmd {

script::script()
  : L(luaL_newstate()) //< create Lua state
{
    // load Lua standard libraries
    luaL_openlibs(L);
    // set Lua package path
    package_path();
    // print Lua stack trace on error
    luabind::set_pcall_callback(&script::traceback);
    // translate C++ standard exceptions into error messages
    register_exception_handlers();
    // load HALMD Lua C++ wrapper
    luaopen_libhalmd(L);
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
 * Load HALMD Lua library
 */
void script::load_library()
{
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
void script::run()
{
    using namespace luabind;

    try {
        slot_function_type slot = resume_function<slot_function_type>(L, "setup");

        // Some C++ modules are only needed during the Lua script stage,
        // e.g. the trajectory reader. To make sure these modules are
        // destructed before running the simulation, invoke the Lua
        // garbage collector now.
        lua_gc(L, LUA_GCCOLLECT, 0);

        slot();

        while (lua_status(L) == LUA_YIELD)
        {
            slot_function_type slot = resume<slot_function_type>(L);

            lua_gc(L, LUA_GCCOLLECT, 0);

            slot();
        }
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
    lua_getglobal(L, "debug");
    lua_pushliteral(L, "traceback");
    lua_rawget(L, -2);
    lua_remove(L, -2);
    lua_call(L, 0, 1);
    lua_concat(L, 3);
    return 1;
}

static void translate_logic_error(lua_State* L, logic_error const& e)
{
    lua_pushliteral(L, "Logic error: ");
    lua_pushstring(L, e.what());
    lua_concat(L, 2);
}

static void translate_domain_error(lua_State* L, domain_error const& e)
{
    lua_pushliteral(L, "Domain error: ");
    lua_pushstring(L, e.what());
    lua_concat(L, 2);
}

static void translate_invalid_argument(lua_State* L, invalid_argument const& e)
{
    lua_pushliteral(L, "Invalid argument: ");
    lua_pushstring(L, e.what());
    lua_concat(L, 2);
}

static void translate_length_error(lua_State* L, length_error const& e)
{
    lua_pushliteral(L, "Length error: ");
    lua_pushstring(L, e.what());
    lua_concat(L, 2);
}

static void translate_out_of_range(lua_State* L, out_of_range const& e)
{
    lua_pushliteral(L, "Out-of-range error: ");
    lua_pushstring(L, e.what());
    lua_concat(L, 2);
}

static void translate_runtime_error(lua_State* L, runtime_error const& e)
{
    lua_pushliteral(L, "Runtime error: ");
    lua_pushstring(L, e.what());
    lua_concat(L, 2);
}

static void translate_range_error(lua_State* L, range_error const& e)
{
    lua_pushliteral(L, "Range error: ");
    lua_pushstring(L, e.what());
    lua_concat(L, 2);
}

static void translate_overflow_error(lua_State* L, overflow_error const& e)
{
    lua_pushliteral(L, "Overflow error: ");
    lua_pushstring(L, e.what());
    lua_concat(L, 2);
}

static void translate_underflow_error(lua_State* L, underflow_error const& e)
{
    lua_pushliteral(L, "Underflow error: ");
    lua_pushstring(L, e.what());
    lua_concat(L, 2);
}

void script::register_exception_handlers()
{
    // C++ standard exceptions
    // http://www.cplusplus.com/reference/std/stdexcept/
    using namespace luabind;
    register_exception_handler<logic_error>(&translate_logic_error);
    register_exception_handler<domain_error>(&translate_domain_error);
    register_exception_handler<invalid_argument>(&translate_invalid_argument);
    register_exception_handler<length_error>(&translate_length_error);
    register_exception_handler<out_of_range>(&translate_out_of_range);
    register_exception_handler<runtime_error>(&translate_runtime_error);
    register_exception_handler<range_error>(&translate_range_error);
    register_exception_handler<overflow_error>(&translate_overflow_error);
    register_exception_handler<underflow_error>(&translate_underflow_error);
}

} // namespace halmd
