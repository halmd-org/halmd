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
#include <luabind/luabind.hpp> // luabind::object
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
    script();
    void dofile(std::string const& file_name);
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
    /** Lua script function */
    luabind::object script_;
};

} // namespace halmd

#endif /* ! HALMD_SCRIPT_HPP */
