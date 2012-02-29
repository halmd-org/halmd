/*
 * Copyright © 2008-2012  Peter Colberg and Felix Höfling
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

#ifndef HALMD_MDSIM_INTEGRATOR_HPP
#define HALMD_MDSIM_INTEGRATOR_HPP

#include <lua.hpp>

namespace halmd {
namespace mdsim {

template <int dimension>
class integrator
{
public:
    static void luaopen(lua_State* L);

    integrator() {}
    virtual ~integrator() {}

    /**
     * update positions and possibly velocities,
     * e.g., first leap-frog step of velocity-Verlet integrator
     */
    virtual void integrate() = 0;

    /**
     * update velocities,
     * e.g., second leap-frog step of velocity-Verlet integrator
     */
    virtual void finalize() = 0;
    /** returns integration timestep */
    virtual double timestep() const = 0;
    /** set integration timestep */
    virtual void timestep(double timestep) = 0;
};

} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_INTEGRATOR_HPP */
