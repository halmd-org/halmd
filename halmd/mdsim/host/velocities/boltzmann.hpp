/*
 * Copyright © 2010 Felix Höfling
 * Copyright © 2008-2012 Peter Colberg
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

#ifndef HALMD_MDSIM_HOST_VELOCITIES_BOLTZMANN_HPP
#define HALMD_MDSIM_HOST_VELOCITIES_BOLTZMANN_HPP

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/random/host/random.hpp>
#include <halmd/utility/profiler.hpp>

#include <lua.hpp>

#include <memory>

namespace halmd {
namespace mdsim {
namespace host {
namespace velocities {

template <int dimension, typename float_type>
class boltzmann
{
public:
    typedef host::particle<dimension, float_type> particle_type;
    typedef random::host::random random_type;

    boltzmann(
        std::shared_ptr<particle_type> particle
      , std::shared_ptr<random_type> random
      , double temperature
      , std::shared_ptr<halmd::logger> logger = std::make_shared<halmd::logger>()
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
     * Set temperature.
     */
    void set_temperature(double temperature);

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    typedef typename particle_type::vector_type vector_type;
    typedef typename particle_type::size_type size_type;
    typedef typename particle_type::velocity_array_type velocity_array_type;
    typedef typename particle_type::mass_array_type mass_array_type;

    /** system state */
    std::shared_ptr<particle_type> particle_;
    /** random number generator */
    std::shared_ptr<random_type> random_;
    /** module logger */
    std::shared_ptr<logger> logger_;
    /** temperature */
    float_type temp_;

    typedef utility::profiler::accumulator_type accumulator_type;
    typedef utility::profiler::scoped_timer_type scoped_timer_type;

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
