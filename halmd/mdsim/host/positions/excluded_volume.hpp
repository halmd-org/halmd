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

#ifndef HALMD_MDSIM_HOST_POSITIONS_EXCLUDED_VOLUME_HPP
#define HALMD_MDSIM_HOST_POSITIONS_EXCLUDED_VOLUME_HPP

#include <boost/make_shared.hpp>
#include <boost/multi_array.hpp>
#include <lua.hpp>
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

    static void luaopen(lua_State* L);

    excluded_volume(
        boost::shared_ptr<box_type const> box
      , float_type cell_length
      , boost::shared_ptr<logger_type> logger = boost::make_shared<logger_type>()
    );
    void exclude_sphere(
        vector_type const& centre
      , float_type diameter
    );
    void exclude_spheres(
        sample_type const& sample
      , std::vector<float_type> diameter
    );
    bool place_sphere(
        vector_type const& centre
      , float_type diameter
    );

private:
    /** simulation box */
    boost::shared_ptr<box_type const> box_;
    /** module logger */
    boost::shared_ptr<logger_type> logger_;
};

} // namespace mdsim
} // namespace host
} // namespace positions
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_POSITIONS_EXCLUDED_VOLUME_HPP */
