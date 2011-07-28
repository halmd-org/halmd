/*
 * Copyright © 2010-2011  Peter Colberg
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

#ifndef HALMD_SCRIPT_HPP
#define HALMD_SCRIPT_HPP

#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>
#include <lua.hpp>
#include <stdexcept>
#include <string>

#include <halmd/utility/options_parser.hpp>
#include <halmd/utility/signal.hpp>

namespace halmd {

/**
 * HALMD scripting engine
 */
class script
  : boost::noncopyable
{
public:
    /**
     * Lua thread
     *
     * This class provides a RAII-safe wrapper around a Lua thread.
     * Lua threads allow functions, or coroutines, to run concurrently.
     * In contrast to operating systems threads, however, only one
     * function is running at any time. While one function is running,
     * all other functions wait for it to yield or return. A function
     * that has yielded continues execution after being resumed by
     * another function.
     *
     * http://www.lua.org/pil/9.html
     * http://www.lua.org/manual/5.1/manual.html#2.11
     */
    class thread
      : boost::noncopyable
    {
    public:
        /**
         * Create a Lua thread with lua_newthread, which pushes the new
         * thread onto the stack. Create a reference to the thread, which
         * also pops the thread from the stack. The reference ensures that
         * the thread is not garbage collected.
         */
        explicit thread(lua_State* L) : L(lua_newthread(L)) , L_(L)
        {
            ref_ = luaL_ref(L_, LUA_REGISTRYINDEX);
        }

        /**
         * Destroy the reference to the thread, which releases the thread
         * for garbage collection.
         */
        ~thread()
        {
            luaL_unref(L_, LUA_REGISTRYINDEX, ref_);
        }

        /**
         * Lua state of the thread. The thread state shares the global
         * variables with the master state, but has its own execution
         * stack.
         */
        lua_State* const L;

    private:
        /** master Lua state */
        lua_State* L_;
        /** reference to Lua thread */
        int ref_;
    };

    script();
    void dofile(std::string const& file_name);
    void load_library();
    void options(options_parser& parser);
    void parsed(po::variables_map const& vm);
    void run();

    static int traceback(lua_State* L);

    //! Lua state
    // Expose Lua state for convenient use in unit tests.
    lua_State* const L;

private:
    typedef signal<void ()>::slot_function_type slot_function_type;

    void package_path();
    static void register_exception_handlers();

    /** RAII wrapper of Lua state */
    boost::shared_ptr<lua_State const> const L_;
};

} // namespace halmd

#endif /* ! HALMD_SCRIPT_HPP */
