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

#ifndef HALMD_MDSIM_FORCE_HPP
#define HALMD_MDSIM_FORCE_HPP

#include <lua.hpp>

namespace halmd {
namespace mdsim {

/**
 * The force module computes all interparticle forces.
 * The current particle positions are read from and
 * the result is stored in the particle module. Periodic
 * boundaries may be taken into account by reference to the
 * box module.
 */

template <int dimension>
class force
{
public:
    static void luaopen(lua_State* L);

    force() {}
    virtual ~force() {}
    virtual void compute() = 0;
    virtual void aux_enable() = 0;
};

} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_FORCE_HPP */
