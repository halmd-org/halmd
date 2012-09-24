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

#ifndef HALMD_MDSIM_GPU_NEIGHBOUR_HPP
#define HALMD_MDSIM_GPU_NEIGHBOUR_HPP

#include <halmd/utility/cache.hpp>

#include <lua.hpp>

#include <cuda_wrapper/cuda_wrapper.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {

/**
 * GPU neighbour lists interface
 *
 * This class provides implementation-independent access to on-GPU
 * neighbour lists for force modules with truncated potentials.
 */
class neighbour
{
public:
    typedef cuda::vector<unsigned int> array_type;

    virtual ~neighbour() {}
    /** Lua bindings */
    static void luaopen(lua_State* L);
    /** neighbour lists */
    virtual cache<array_type> const& g_neighbour() = 0;
    /** number of placeholders per neighbour list */
    virtual unsigned int size() const = 0;
    /** neighbour list stride */
    virtual unsigned int stride() const = 0;
};

} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_NEIGHBOUR_HPP */
