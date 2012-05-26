/*
 * Copyright © 2010 Felix Höfling
 * Copyright © 2008-2012 Peter Colberg
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

#ifndef HALMD_MDSIM_HOST_VELOCITIES_BOLTZMANN_HPP
#define HALMD_MDSIM_HOST_VELOCITIES_BOLTZMANN_HPP

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <lua.hpp>
#include <utility>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/host/velocity.hpp>
#include <halmd/random/host/random.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd {
namespace mdsim {
namespace host {
namespace velocities {

template <int dimension, typename float_type>
class boltzmann
  : public host::velocity<dimension, float_type>
{
public:
    typedef host::particle<dimension, float_type> particle_type;
    typedef random::host::random random_type;
    typedef logger logger_type;

    boltzmann(
        boost::shared_ptr<particle_type> particle
      , boost::shared_ptr<random_type> random
      , double temperature
      , boost::shared_ptr<logger_type> logger = boost::make_shared<logger_type>()
    );

    /**
     * Initialise velocities from Maxwell-Boltzmann distribution
     */
    void set();

    /**
     * Returns temperature.
     */
    float_type temperature() const
    {
        return temp_;
    }

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    typedef host::velocity<dimension, float_type> _Base;
    typedef typename particle_type::vector_type vector_type;

    /**
     * Assign new velocities from Gaussian distribution
     *
     * @param sigma width of distribution
     * @returns mean velocity and mean-square velocity
     */
    std::pair<vector_type, float_type> gaussian(float_type sigma);

    /** system state */
    boost::shared_ptr<particle_type> particle_;
    /** random number generator */
    boost::shared_ptr<random_type> random_;
    /** module logger */
    boost::shared_ptr<logger_type> logger_;
    /** temperature */
    float_type temp_;

    typedef utility::profiler profiler_type;
    typedef typename profiler_type::accumulator_type accumulator_type;
    typedef typename profiler_type::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type set;
    };

    /** profiling runtime accumulators */
    runtime runtime_;
};

} // namespace velocities
} // namespace host
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_VELOCITIES_BOLTZMANN_HPP */
