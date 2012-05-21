/*
 * Copyright © 2008-2012  Peter Colberg and Felix Höfling
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

#ifndef HALMD_MDSIM_HOST_INTEGRATORS_VERLET_NVT_ANDERSEN_HPP
#define HALMD_MDSIM_HOST_INTEGRATORS_VERLET_NVT_ANDERSEN_HPP

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <lua.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/random/host/random.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd {
namespace mdsim {
namespace host {
namespace integrators {

template <int dimension, typename float_type>
class verlet_nvt_andersen
{
public:
    typedef host::particle<dimension, float_type> particle_type;
    typedef mdsim::box<dimension> box_type;
    typedef random::host::random random_type;
    typedef logger logger_type;

private:
    typedef typename particle_type::vector_type vector_type;

public:
    /**
     * Initialise Verlet-Andersen integrator.
     */
    verlet_nvt_andersen(
        boost::shared_ptr<particle_type> particle
      , boost::shared_ptr<box_type const> box
      , boost::shared_ptr<random_type> random
      , float_type timestep
      , float_type temperature
      , float_type coll_rate
      , boost::shared_ptr<logger_type> logger = boost::make_shared<logger_type>()
    );

    /**
     * First leapfrog half-step of velocity-Verlet algorithm
     */
    void integrate();

    /**
     * Second leapfrog half-step of velocity-Verlet algorithm
     */
    void finalize();

    /**
     * Set integration time-step.
     */
    void set_timestep(double timestep);

    /**
     * Returns integration time-step.
     */
    double timestep() const
    {
        return timestep_;
    }

    /**
     * Set temperature of heat bath.
     */
    void set_temperature(double temperature);

    /**
     * Returns temperature of heat bath.
     */
    double temperature() const
    {
        return temperature_;
    }

    /**
     * Returns collision rate with the heat bath.
     */
    float_type collision_rate() const
    {
        return coll_rate_;
    }

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    /** system state */
    boost::shared_ptr<particle_type> particle_;
    /** simulation domain */
    boost::shared_ptr<box_type const> box_;
    /** random number generator */
    boost::shared_ptr<random_type> random_;
    /** integration time-step */
    float_type timestep_;
    /** half time-step */
    float_type timestep_half_;
    /** temperature of the heat bath */
    float_type temperature_;
    /** square root of temperature */
    float_type sqrt_temperature_;
    /** collision rate with the heat bath */
    float_type coll_rate_;
    /** probability of a collision with the heat bath during a timestep */
    float_type coll_prob_;
    /** module logger */
    boost::shared_ptr<logger_type> logger_;

    typedef utility::profiler profiler_type;
    typedef typename profiler_type::accumulator_type accumulator_type;
    typedef typename profiler_type::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type integrate;
        accumulator_type finalize;
    };

    /** profiling runtime accumulators */
    runtime runtime_;
};

} // namespace integrators
} // namespace host
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_INTEGRATORS_VERLET_NVT_ANDERSEN_HPP */
