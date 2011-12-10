/*
 * Copyright © 2011  Peter Colberg
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

#ifndef HALMD_MDSIM_HOST_NEIGHBOUR_HPP
#define HALMD_MDSIM_HOST_NEIGHBOUR_HPP

#include <lua.hpp>
#include <vector>

#include <halmd/mdsim/neighbour.hpp>

namespace halmd {
namespace mdsim {
namespace host {

/**
 * host neighbour lists interface
 *
 * This class provides implementation-independent access to host
 * neighbour lists for force modules with truncated potentials.
 */
class neighbour
  : public mdsim::neighbour
{
public:
    typedef std::vector<unsigned int> neighbour_list;

    /** Lua bindings */
    static void luaopen(lua_State* L);
    /** neighbour lists */
    virtual std::vector<neighbour_list> const& lists() const = 0;
};

} // namespace host
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_NEIGHBOUR_HPP */
