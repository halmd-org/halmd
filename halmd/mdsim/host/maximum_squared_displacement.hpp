/*
 * Copyright © 2008-2011  Peter Colberg
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

#ifndef HALMD_MDSIM_HOST_MAXIMUM_SQUARED_DISPLACEMENT_HPP
#define HALMD_MDSIM_HOST_MAXIMUM_SQUARED_DISPLACEMENT_HPP

#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/shared_ptr.hpp>
#include <lua.hpp>
#include <vector>

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/force.hpp>
#include <halmd/mdsim/host/particle.hpp>

namespace halmd {
namespace mdsim {
namespace host {

template <int dimension, typename float_type>
class maximum_squared_displacement
{
public:
    typedef host::particle<dimension, float_type> particle_type;
    typedef typename particle_type::vector_type vector_type;
    typedef mdsim::box<dimension> box_type;


    static void luaopen(lua_State* L);

    maximum_squared_displacement(
        boost::shared_ptr<particle_type const> particle
      , boost::shared_ptr<box_type const> box
    );
    void zero();
    float_type compute();

private:
    //! system state
    boost::shared_ptr<particle_type const> particle_;
    //! simulation box
    boost::shared_ptr<box_type const> box_;
    /* particle positions at last maximum_squared_displacement list update */
    std::vector<vector_type> r0_;
};

} // namespace host
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_MAXIMUM_SQUARED_DISPLACEMENT_HPP */
