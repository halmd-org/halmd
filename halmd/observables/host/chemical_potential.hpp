/*
 * Copyright © 2016 Felix Höfling
 * Copyright © 2016 Arthur Straube
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

#ifndef HALMD_OBSERVABLES_HOST_CHEMICAL_POTENTIAL_HPP
#define HALMD_OBSERVABLES_HOST_CHEMICAL_POTENTIAL_HPP

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/mdsim/host/positions/lattice.hpp> // FIXME random.hpp
#include <halmd/numeric/accumulator.hpp>
#include <halmd/utility/cache.hpp>
#include <halmd/utility/profiler.hpp>

#include <lua.hpp>
#include <memory>
#include <vector>

namespace halmd {
namespace observables {
namespace host {

/**
 * Compute the excess chemical potential using Widom's insertion method.
 */

template <int dimension, typename float_type>
class chemical_potential
{
public:
    typedef mdsim::host::particle<dimension, float_type> particle_type;
    typedef mdsim::host::positions::lattice<dimension, float_type> position_type;
    typedef std::vector<halmd::accumulator<double>> result_type;

    typedef typename particle_type::vector_type vector_type;
    typedef typename particle_type::size_type size_type;

    static void luaopen(lua_State* L);

    chemical_potential(
        std::shared_ptr<particle_type> test_particle
      , std::shared_ptr<position_type> position
      , double temperature
      , std::vector<unsigned int> ntest_particle
      , std::shared_ptr<halmd::logger> logger = std::make_shared<halmd::logger>()
    );

    /**
     * Perform Widom's insertion method and compute excess chemical potential.
     *
     * Return excess chemical potential per species as triple (mean, error of mean, count).
     */
    result_type const& sample();

    /**
     * Re-assign random positions to test particles.
     */
    void set_position()
    {
        position_->set();
    }

    /**
     * Set temperature.
     */
    void set_temperature(double temperature);

    /**
     * Returns temperature.
     */
    double temperature() const
    {
        return temperature_;
    }

    /**
     * Returns test particle instance.
     */
    std::shared_ptr<particle_type const> test_particle() const
    {
        return test_particle_;
    }

    /**
     * length of result vector, i.e., the number of particles species'
     */
    unsigned int result_size() const
    {
        return mu_ex_.size();
    }

private:
    typedef typename particle_type::en_pot_array_type en_pot_array_type;

    /** test particles */
    std::shared_ptr<particle_type> test_particle_;
    /** module to set particle positions */
    std::shared_ptr<position_type> position_;

    /** number of test particles per species */
    std::vector<unsigned int> ntest_particle_;
    /** temperature */
    double temperature_;

    /** excess chemical potential per species */
    result_type mu_ex_;
    /** cache observer of excess chemical potential */
    cache<> mu_ex_cache_;

    /** module logger */
    std::shared_ptr<logger> logger_;

    typedef halmd::utility::profiler::accumulator_type accumulator_type;
    typedef halmd::utility::profiler::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type sample;
    };

    /** profiling runtime accumulators */
    runtime runtime_;
};

} // namespace host
} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_CHEMICAL_POTENTIAL_HPP */
