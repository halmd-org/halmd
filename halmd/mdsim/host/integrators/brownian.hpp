/*
 * Copyright ©  2015 Manuel Dibak
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

#ifndef HALMD_MDSIM_HOST_INTEGRATORS_BROWNIAN_HPP
#define HALMD_MDSIM_HOST_INTEGRATORS_BROWNIAN_HPP

#include <lua.hpp>
#include <memory>

#include <boost/numeric/ublas/vector.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/host/particle.hpp>
#include <halmd/utility/profiler.hpp>
#include <halmd/random/host/random.hpp>

namespace halmd {
namespace mdsim {
namespace host {
namespace integrators {

template <int dimension, typename float_type>
class brownian
{
public:
    typedef host::particle<dimension, float_type> particle_type;
    typedef mdsim::box<dimension> box_type;
    typedef random::host::random random_type;

    typedef typename particle_type::vector_type vector_type;
    typedef boost::numeric::ublas::matrix<float_type> matrix_type;

    brownian(
        std::shared_ptr<particle_type> particle
      , std::shared_ptr<random_type> random
      , std::shared_ptr<box_type const> box
      , float_type timestep
      , float_type temperature
      , const float_type& diff_const
      , std::shared_ptr<halmd::logger> logger = std::make_shared<halmd::logger>()
    );

    /**
     * Brownian integration step.
     */
    void integrate();

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
     * Bind class to Lua.
     */
static void luaopen(lua_State* L);

private:
    typedef typename particle_type::position_array_type position_array_type;
    typedef typename particle_type::image_array_type image_array_type;
    typedef typename particle_type::velocity_array_type velocity_array_type;
    typedef typename particle_type::size_type size_type;

    typedef utility::profiler::accumulator_type accumulator_type;
    typedef utility::profiler::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type integrate;
    };

    void update_displacement(
      float_type diff_const
      , vector_type& r
      , vector_type const& f
      , vector_type const& eta
    );
    
    /**
     * set rng and random seed
     */
    unsigned int seed_;

    std::shared_ptr<particle_type> particle_;
    std::shared_ptr<random_type> random_;
    std::shared_ptr<box_type const> box_;

    /** integration time-step */
    float_type timestep_;
    /** temperature */
    float_type temperature_;
    /** profiling runtime accumulators */
    runtime runtime_;
    /** diffusion constant */
    float_type diff_const_;
    /** module logger */
    std::shared_ptr<logger> logger_;
};

} // namespace integrators
} // namespace host
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_HOST_INTEGRATORS_BROWNIAN_HPP */
