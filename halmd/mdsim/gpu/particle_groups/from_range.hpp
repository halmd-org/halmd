/*
 * Copyright © 2012 Peter Colberg
 * Copyright © 2012 Felix Höfling
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

#ifndef HALMD_MDSIM_GPU_PARTICLE_GROUPS_FROM_RANGE_HPP
#define HALMD_MDSIM_GPU_PARTICLE_GROUPS_FROM_RANGE_HPP

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
 * Select particles of a given particle instance by a contiguous range of particle tags.
 */
template <typename particle_type>
class from_range
  : public particle_group
{
public:
    typedef typename particle_group::array_type array_type;
    typedef typename particle_group::size_type size_type;
    typedef std::pair<size_type, size_type> range_type;
    typedef logger logger_type;

    /**
     * Select by tag range [begin, end).
     */
    from_range(
        std::shared_ptr<particle_type const> particle
      , range_type const& range
      , std::shared_ptr<logger_type> logger = std::make_shared<logger_type>()
    );

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
    /** validate tag range */
    range_type const& check_range(range_type const&);

    /** particle instance */
    std::shared_ptr<particle_type const> const particle_;
    /** particle tag range */
    range_type const range_;
    /** module logger */
    std::shared_ptr<logger_type> logger_;
    /** ordered sequence of particle indices */
    cache<array_type> ordered_;
    /** unordered sequence of particle indices */
    cache<array_type> unordered_;
    /** number of particles */
    cache<size_type> size_;
    /** cache observer of particle reverse tags */
    cache<> ordered_cache_;
    /** cache observer of particle reverse tags */
    cache<> unordered_cache_;
};

} // namespace particle_groups
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_PARTICLE_GROUPS_FROM_RANGE_HPP */
