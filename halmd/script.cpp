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

#include <luabind/class_info.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/script.hpp>
#include <halmd/utility/filesystem.hpp>
#include <halmd/utility/lua/lua.hpp>
#include <halmd/utility/lua/program_options.hpp>
#include <halmd/version.h>

using namespace boost;
using namespace std;

HALMD_LUA_API int luaopen_halmd_base(lua_State* L);

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
    // set Lua package path
    package_path();
    // set Lua C package path
    package_cpath();
    // print Lua stack trace on error
    luabind::set_pcall_callback(&script::traceback);
    // translate C++ standard exceptions into error messages
    register_exception_handlers();
    // setup global structures and Lua class support
    luabind::open(L);
    // class_info(), class_names()
    luabind::bind_class_info(L);
    // load HALMD Lua C++ wrapper
    luaopen_halmd_base(L);
    // prepare Lua 5.2 compatible environment
    lua_compat();
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
    using namespace luabind;

#if LUA_VERSION_NUM < 502
    // function unpack was moved into the table library
    globals(L)["table"]["unpack"] = globals(L)["unpack"];
    globals(L)["unpack"] = nil;
#endif
}

/**
 * Load HALMD Lua library
 */
void script::load_library()
{
    using namespace luabind;

    try {
        script_ = call_function<object>(L, "require", "halmd.default");
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
    using namespace luabind;

    // error handler passed to lua_pcall as last argument
    lua_pushcfunction(L, &script::traceback);

    if (luaL_loadfile(L, file_name.c_str()) || lua_pcall(L, 0, 1, 1)) {
        LOG_ERROR(lua_tostring(L, -1));
        lua_pop(L, 1); //< remove error message
        throw runtime_error("failed to load Lua script");
    }

    // store return value of script as HALMD script function
    // this function will be called later in script::run,
    // after the command-line options have been parsed
    script_ = object(from_stack(L, -1));
    lua_pop(L, 1);
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
 *
 * While the simulation is setup with Lua scripting, it entirely runs in
 * C++, outside of the Lua interpreter. This avoids unnecessary Lua errors
 * in case a C++ exception is thrown inside a module, e.g. if a defect GPU
 * fails during a long simulation run and a CUDA error is thrown. Further
 * it simplifies stack tracebacks when running inside a debugger.
 *
 * To allow the script to run multiple partial simulations, we use Lua
 * coroutines. The Lua function script() is started, and executes until
 * yielding a slot, e.g. sampler.setup or sampler.run. The slot is then
 * executed outside of the Lua interpreter. The Lua function is then
 * resumed. Script execution finishes when the Lua function returns.
 *
 * We use the Lua C API instead of luabind::resume_function and
 * luabind::resume, as the latter do not handle errors properly,
 * but abort with an assertion error (or segmentation fault if
 * compiled with NDEBUG).
 */
void script::run()
{
    using namespace luabind;

    // create a new Lua thread
    //
    // Lua threads allow functions, or coroutines, to run concurrently.
    // In contrast to operating systems threads, however, only one
    // function is running at any time. While one function is running,
    // all other functions wait for it to yield or return. A function
    // that has yielded continues execution after being resumed by
    // another function.
    //
    // http://www.lua.org/pil/9.html
    // http://www.lua.org/manual/5.1/manual.html#2.11
    //
    // override member pointer L to master state with thread state,
    // to prevent errors due to accidental use of master state
    lua_State* const L = lua_newthread(script::L);

    // The object wrapper ensures that the thread is released for garbage
    // collection in case a C++ exception is thrown. Note that the new
    // thread is immediately popped from the stack.
    object thread(from_stack(script::L, -1));
    lua_pop(script::L, 1);

    // push HALMD script function onto stack
    script_.push(L);
    if (lua_isnil(L, -1)) {
        throw runtime_error("missing callable return value from HALMD script");
    }

    int status;
    do {
        // if Lua function is on top of the stack, create a new coroutine
        // from it, otherwise resume execution of the existing coroutine
#if LUA_VERSION_NUM < 502
        status = lua_resume(L, 0);
#else
        status = lua_resume(L, NULL, 0);
#endif

        // lua_resume returns
        //  - 0 if the function has returned successfully
        //  - LUA_YIELD if the function has yielded
        //  - other values if an error occurred
        if (status == 0) {
            // we expect no return value
        }
        else if (status == LUA_YIELD) {
            // we expect a slot as the yield value
            object ret(from_stack(L, -1));
            lua_pop(L, 1);
            slot_function_type slot = object_cast<slot_function_type>(ret);

            // Some C++ modules are only needed during the Lua script stage,
            // e.g. the trajectory reader. To make sure these modules are
            // destructed before running the simulation, invoke the Lua
            // garbage collector now.
            lua_gc(L, LUA_GCCOLLECT, 0);

            // execute the slot, e.g. sampler.setup or sampler.run
            slot();
        }
        else {
            script::traceback(L);
            LOG_ERROR(lua_tostring(L, -1));
            lua_pop(L, 1); // remove error message
            throw runtime_error("failed to run simulation script");
        }
    }
    while (status == LUA_YIELD);
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
