/*
 * Copyright © 2011 Michael Kopp
 * Copyright © 2014 Felix Höfling
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

#ifndef HALMD_MDSIM_HOST_MOBILITIES_OSEEN_HPP
#define HALMD_MDSIM_HOST_MOBILITIES_OSEEN_HPP

#include <boost/shared_ptr.hpp>
#include <lua.hpp>

#include <halmd/mdsim/box.hpp> // reduce_periodic
#include <halmd/mdsim/mobility.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd {
namespace mdsim {
namespace host {
namespace mobilities {

template <int dimension, typename float_type>
class oseen
  : public mdsim::mobility<dimension>
{
public:
    typedef mdsim::mobility<dimension> _Base;
    typedef host::particle<dimension, float_type> particle_type;
    typedef typename particle_type::vector_type vector_type;
    typedef typename mdsim::box<dimension> box_type;

    static void luaopen(lua_State* L);

    /**
     * @param particle shared pointer to particle module describing current microscopic state
     * @param box the particles reside in to evaluate periodic boundary conditions
     * @param radius hydrodymanic radius of particle (global rad. for each particle)
     * @param viscosity dynamic visosity of fluid
     * @param order of accuracy of hydrodynamic interactions in (a/r) (1,2: Oseen, 3: Rotne Prager)
     */
    oseen(
        boost::shared_ptr<particle_type> particle
      , boost::shared_ptr<box_type> box
      , float radius
      , float viscosity
      , int order
      , boost::shared_ptr<halmd::logger> logger = boost::make_shared<halmd::logger>()
    );

    // inherited functions
    virtual void compute_velocity();
    virtual void compute();

    /** returns hydrodynamic radius */
    inline float radius() const
    {
        return radius_;
    }

    /** returns dynamic viscosity */
    inline float viscosity() const
    {
        return viscosity_;
    }

    /** return order of accuracy of hydrodynamic interaction in (a/r) */
    inline int order() const
    {
        return order_;
    }

private:
    /** particle instance */
    boost::shared_ptr<particle_type> particle_;
    /** box instance */
    boost::shared_ptr<box_type> box_;
    /** module logger */
    boost::shared_ptr<logger> logger_;

    /** hydrodynamic radius */
    float radius_;
    /** dynamic viscosity of fluid */
    float viscosity_;
    /** order of accuracy of hydrodynamic interaction in (a/r) */
    int order_;

    struct runtime
    {
        utility::profiler::accumulator_type compute_velocity;
        utility::profiler::accumulator_type compute;
    };

    /** profiling runtime accumulators */
    runtime runtime_;
};

} // namespace mobilities
} // namespace host
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_MOBILITIES_OSEEN_HPP */
