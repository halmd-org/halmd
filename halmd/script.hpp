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

#ifndef HALMD_SCRIPT_HPP
#define HALMD_SCRIPT_HPP

#include <boost/noncopyable.hpp>
#include <lua.hpp>
#include <string>

#include <halmd/utility/options_parser.hpp>

namespace halmd
{

/**
 * HALMD scripting engine
 */
class script
  : boost::noncopyable
{
public:
    script();
    virtual ~script();
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
    void package_path();
    void load_wrapper();
};

} // namespace halmd

#endif /* ! HALMD_SCRIPT_HPP */
