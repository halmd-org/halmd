/*
 * Copyright © 2011  Peter Colberg
 * Copyright © 2012  Nicolas Höft
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

#ifndef HALMD_MDSIM_HOST_POSITIONS_EXCLUDED_VOLUME_HPP
#define HALMD_MDSIM_HOST_POSITIONS_EXCLUDED_VOLUME_HPP

#include <boost/multi_array.hpp>
#include <lua.hpp>
#include <memory>
#include <vector>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/observables/host/samples/phase_space.hpp>

namespace halmd {
namespace mdsim {
namespace host {
namespace positions {

/**
 * Track excluded volume for particle placement
 */
template <int dimension, typename float_type>
class excluded_volume
{
public:
    typedef mdsim::box<dimension> box_type;
    typedef observables::host::samples::phase_space<dimension, float_type> sample_type;
    typedef typename sample_type::vector_type vector_type;
    typedef logger logger_type;

    excluded_volume(
        std::shared_ptr<box_type const> box
      , float_type cell_length
      , std::shared_ptr<logger_type> logger = std::make_shared<logger_type>()
    );

    /**
     * Place a sphere and consider it for future calls to place_sphere()
     */
    void exclude_sphere(
        vector_type const& centre
      , float_type diameter
    );

    /**
     * Place a number of spheres at specified positions with given radii
     */
    void exclude_spheres(
        sample_type const& sample
      , std::vector<float_type> diameter
    );

    /**
     * Test if the placement of a sphere is possible with no overlap with other
     * previously placed spheres
     * 
     * @param center center position of sphere
     * @param diameter diameter of the sphere
     * @return true if the placement is possible
     */
    bool place_sphere(
        vector_type const& centre
      , float_type diameter
    ) const;

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    typedef std::pair<vector_type, float_type> sphere_type;
    typedef fixed_vector<unsigned int, dimension> index_type;
    typedef std::pair<index_type, index_type> index_pair_type;
    typedef std::vector<sphere_type> cell_type;

    index_pair_type sphere_extents(
        vector_type const& centre
      , float_type diameter
    ) const;
    void exclude_sphere_from_cell(
        vector_type const& centre
      , float_type diameter
      , index_type const& index
    );
    bool place_cell(
        vector_type const& centre
      , float_type diameter
      , index_type const& index
    ) const;

    /** simulation box */
    std::shared_ptr<box_type const> box_;
    /** module logger */
    std::shared_ptr<logger_type> logger_;
    /** number of cells per dimension */
    index_type ncell_;
    /** cell edge length per dimension */
    vector_type cell_length_;
    /** n-dimensional cells with excluded spheres */
    boost::multi_array<cell_type, dimension> cell_;
};

} // namespace mdsim
} // namespace host
} // namespace positions
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_POSITIONS_EXCLUDED_VOLUME_HPP */
