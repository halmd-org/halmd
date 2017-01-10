/*
 * Copyright © 2012 Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#ifndef HALMD_MDSIM_GPU_PARTICLE_GROUPS_ALL_HPP
#define HALMD_MDSIM_GPU_PARTICLE_GROUPS_ALL_HPP

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/gpu/particle_group.hpp>
#include <halmd/utility/raw_array.hpp>

#include <lua.hpp>

#include <memory>
#include <utility>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace particle_groups {

/**
 * Select particles of a given particle instance by a contiguous range of particle IDs.
 */
template <typename particle_type>
class all
  : public particle_group
{
public:
    typedef typename particle_group::array_type array_type;
    typedef typename particle_group::size_type size_type;

    /**
     * Select all particles.
     */
    all(std::shared_ptr<particle_type const> particle);

    /**
     * Returns ordered sequence of particle indices.
     */
    virtual cache<array_type> const& ordered();

    /**
     * Returns unordered sequence of particle indices.
     */
    virtual cache<array_type> const& unordered();

    /**
     * Returns number of particles.
     */
    virtual cache<size_type> const& size();

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    /** particle instance */
    std::shared_ptr<particle_type const> const particle_;
    /** unordered sequence of particle indices */
    cache<array_type> unordered_;
    /** number of particles */
    cache<size_type> size_;
};

} // namespace particle_groups
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_PARTICLE_GROUPS_ALL_HPP */
