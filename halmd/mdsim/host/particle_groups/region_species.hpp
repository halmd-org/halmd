/*
 * Copyright © 2014-2015 Nicolas Höft
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

#ifndef HALMD_MDSIM_HOST_PARTICLE_GROUPS_REGION_SPECIES_HPP
#define HALMD_MDSIM_HOST_PARTICLE_GROUPS_REGION_SPECIES_HPP

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/host/particle_group.hpp>
#include <halmd/mdsim/host/particle.hpp>

#include <lua.hpp>
#include <memory>

namespace halmd {
namespace mdsim {
namespace host {
namespace particle_groups {

/**
 * Select particles of a given particle instance from a region in space and of
 * a given species.
 */
template <int dimension, typename float_type, typename geometry_type>
class region_species
  : public particle_group
{
public:
    typedef typename particle_group::array_type array_type;
    typedef typename particle_group::size_type size_type;
    typedef host::particle<dimension, float_type> particle_type;
    typedef typename particle_type::vector_type vector_type;

    enum geometry_selection {
        excluded = 1
      , included = 2
    };

    /**
     * Select by region and species
     */
    region_species(
        std::shared_ptr<particle_type const> particle
      , std::shared_ptr<geometry_type const> geometry
      , geometry_selection geometry_sel
      , unsigned int species
      , std::shared_ptr<halmd::logger> logger = std::make_shared<halmd::logger>()
    );

    /**
     * Returns indices of particles that are within the defined region of the
     * simulation box and of the given species
     */
    cache<array_type> const& selection();

    /**
     * Returns boolean mask for each particle whether it is within the defined
     * region of the simulation box and of the given species
     */
    cache<array_type> const& mask();

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
    typedef typename particle_type::position_array_type position_array_type;

    void update_();

    /** particle instance */
    std::shared_ptr<particle_type const> particle_;
    /** region the particles are sorted by */
    std::shared_ptr<geometry_type const> geometry_;
    geometry_selection geometry_selection_;
    /** select particles of given species */
    unsigned int species_;
    /** module logger */
    std::shared_ptr<logger> logger_;

    /**
     * mask for particles that determines whether they are in-/outside the
     * region, each element has the value 0 or 1, where 0 means outside the
     * region and 1 included in region. This mask is ordered by particle id.
     */
    cache<array_type> mask_;
    /** particle indices of particles in the region and of the species */
    cache<array_type> selection_;
    /** cache observer of position updates for mask */
    cache<> mask_cache_;

    /** ordered sequence of particle indices */
    cache<array_type> ordered_;
    /** number of particles */
    cache<size_type> size_;
    /** cache observer of selection */
    cache<> ordered_cache_;
    /** cache observer of size */
    cache<> size_cache_;

    typedef utility::profiler::scoped_timer_type scoped_timer_type;
    typedef utility::profiler::accumulator_type accumulator_type;

    struct runtime
    {
        accumulator_type update_mask;
        accumulator_type update_selection;
        accumulator_type sort_selection;
    };

    /** profiling runtime accumulators */
    runtime runtime_;
};

} // namespace particle_groups
} // namespace host
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_PARTICLE_GROUPS_REGION_SPECIES_HPP */
