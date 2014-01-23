/*
 * Copyright © 2011 Michael Kopp
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

#ifndef HALMD_MDSIM_MOBILITY_HPP
#define HALMD_MDSIM_MOBILITY_HPP

#include <lua.hpp>

namespace halmd {
namespace mdsim {

/**
 * The mobility module will one day compute the mobility
 * matrix via compute(). Velocities can then be obtained
 * from the matrix product with the forces.
 *  NOT IMPLEMENTED YET!
 * The function compute_velocities() computes the velocities
 * directly and stores them in the particle module. The
 * positions necessary to compute the mobility are read from
 * the particle module.
 */

template <int dimension>
class mobility
{
public:
    static void luaopen(lua_State* L);

    mobility() {}
    virtual ~mobility() {}
    //! Compute mobility matrix. (Not yet implemented)
    virtual void compute() = 0; // matrix
    //! Compute velocities directly from positions and forces.
    virtual void compute_velocities() = 0;
};

} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_MOBILITY_HPP */
