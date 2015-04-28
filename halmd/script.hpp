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

#ifndef HALMD_SCRIPT_HPP
#define HALMD_SCRIPT_HPP

#include <boost/noncopyable.hpp>
#include <lua.hpp>
#include <memory>
#include <string>

namespace halmd {

/**
 * HALMD scripting engine
 */
class script
  : boost::noncopyable
{
public:
    script();
    void dofile(std::string const& filename = std::string());

    static int traceback(lua_State* L);

    //! Lua state
    // Expose Lua state for convenient use in unit tests.
    lua_State* const L;

private:
    void load_luaponte();
    void prepend_package_path(std::string const& path);
    void prepend_package_cpath(std::string const& path);

    /** RAII wrapper of Lua state */
    std::shared_ptr<lua_State const> const L_;
};

} // namespace halmd

#endif /* ! HALMD_SCRIPT_HPP */
