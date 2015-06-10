/*
 * Copyright © 2014-2015 Nicolas Höft
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

#ifndef HALMD_MDSIM_HOST_REGION_HPP
#define HALMD_MDSIM_HOST_REGION_HPP

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/utility/profiler.hpp>

#include <lua.hpp>

#include <memory>
#include <vector>

namespace halmd {
namespace mdsim {
namespace host {

class region_base
{
public:
    typedef std::vector<unsigned int> array_type;
    typedef typename array_type::value_type size_type;

    /**
     * Returns list of particle indices that are in the
     * defined region of the system
     */
    virtual cache<array_type> const& selection() = 0;

    /**
     * Number of particles in the region
     */
    virtual size_type size() = 0;

    /**
     * mask of that specifies if a particle is within the region
     * or outside. Sorted by particle id
     */
    virtual cache<array_type> const& mask() = 0;

    /**
     * Bind class to Lua
     */
    static void luaopen(lua_State* L);
};

template <int dimension, typename float_type, typename geometry_type>
class region
  : public region_base
{
public:
    typedef region_base::array_type array_type;
    typedef region_base::size_type size_type;
    typedef host::particle<dimension, float_type> particle_type;
    typedef typename particle_type::vector_type vector_type;
    typedef mdsim::box<dimension> box_type;

    enum geometry_selection {
        excluded = 1
      , included = 2
    };

    /**
     * Bind class to Lua
     */
    static void luaopen(lua_State* L);

    region(
        std::shared_ptr<particle_type const> particle
      , std::shared_ptr<box_type const> box
      , std::shared_ptr<geometry_type const> geometry
      , geometry_selection geometry_sel
      , std::shared_ptr<halmd::logger> logger = std::make_shared<halmd::logger>()
    );

    /**
     * Returns particle indices that are within the
     * defined region of the simulation box
     */
    cache<array_type> const& selection();

    size_type size()
    {
        update_();
        return selection_->size();
    }

    cache<array_type> const& mask();

private:
    typedef typename particle_type::position_array_type position_array_type;
    typedef typename particle_type::position_type position_type;

    void update_();

    //! system state
    std::shared_ptr<particle_type const> particle_;
    //! simulation box
    std::shared_ptr<box_type const> box_;
    /** module logger */
    std::shared_ptr<logger> logger_;
    //! region the particles are sorted by
    std::shared_ptr<geometry_type const> geometry_;
    /** cache observer of position updates for mask */
    cache<> mask_cache_;

    geometry_selection geometry_selection_;
    /**
     * mask for particles that determines whether they are in-/outside the region,
     * each element has the value 0 or 1, where 0 means outside
     * the region and 1 included in region.
     * This mask is ordered by particle id.
     */
    cache<array_type> mask_;
    /** particle indices of particles in the region */
    cache<array_type> selection_;

    typedef utility::profiler::scoped_timer_type scoped_timer_type;
    typedef utility::profiler::accumulator_type accumulator_type;

    struct runtime
    {
        accumulator_type update_mask;
        accumulator_type update_selection;
    };
    /** profiling runtime accumulators */
    runtime runtime_;
};

/**
 * Copy particle ids of particles in region into given array
 */
template <typename region_type, typename iterator_type>
inline iterator_type
get_selection(region_type& region, iterator_type const& first)
{
    auto const& selection = read_cache(region.selection());
    return std::copy(selection.begin(), selection.end(), first);
}

/**
 * Copy region mask into given array
 */
template <typename region_type, typename iterator_type>
inline iterator_type
get_mask(region_type& region, iterator_type const& first)
{
    auto const& mask = read_cache(region.mask());
    return std::copy(mask.begin(), mask.end(), first);
}

} // namespace host
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_REGION_HPP */
