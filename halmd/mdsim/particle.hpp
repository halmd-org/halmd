/*
 * Copyright © 2008-2010  Peter Colberg
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

#ifndef HALMD_MDSIM_PARTICLE_HPP
#define HALMD_MDSIM_PARTICLE_HPP

#include <lua.hpp>
#include <vector>

#include <halmd/options.hpp>

namespace halmd
{
namespace mdsim
{

/**
 * The particle module describes the current microsopic state of the system.
 * It holds positions, velocities, forces, particle identities
 * ... and the neighbour lists.
 */

template <int dimension>
class particle
{
public:
    static void options(po::options_description& desc);

    static void luaopen(lua_State* L);

    particle(std::vector<unsigned int> const& particles);
    virtual ~particle() {}
    virtual void set() = 0;

    /** number of particles in simulation box */
    unsigned int const nbox;
    /** number of particle types */
    unsigned int const ntype;
    /** number of particles per type */
    std::vector<unsigned int> const ntypes;
};

} // namespace mdsim

} // namespace halmd

#endif /* ! HALMD_MDSIM_PARTICLE_HPP */
