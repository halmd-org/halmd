/*
 * Copyright © 2015 Manuel Dibak
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

#ifndef HALMD_MDSIM_GPU_INTEGRATORS_BROWNIAN_HPP
#define HALMD_MDSIM_GPU_INTEGRATORS_BROWNIAN_HPP

#include <lua.hpp>
#include <memory>

#include <boost/numeric/ublas/vector.hpp>
#include <cuda_wrapper/cuda_wrapper.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/gpu/integrators/brownian_kernel.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/random/gpu/random.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace integrators {

template <int dimension, typename float_type, typename RandomNumberGenerator>
class brownian
{
public:
    typedef gpu::particle<dimension, float_type> particle_type;
    typedef mdsim::box<dimension> box_type;
    typedef typename particle_type::vector_type vector_type;
    typedef boost::numeric::ublas::matrix<float_type> matrix_type;
    typedef random::gpu::random<RandomNumberGenerator> random_type;
    typedef typename random_type::rng_type rng_type;
    typedef brownian_wrapper<dimension, float_type, rng_type> wrapper_type;

    static void luaopen(lua_State* L);

    brownian(
        std::shared_ptr<particle_type> particle
      , std::shared_ptr<random_type> random
      , std::shared_ptr<box_type const> box
      , double timestep
      , double temperature
      , matrix_type const& diff_const
      , std::shared_ptr<halmd::logger> logger = std::make_shared<halmd::logger>()
    );

    void integrate();

    /**
     * Set integration time-step.
     */
    void set_timestep(double timestep);

    /**
     * Set integration time-step.
     */
void set_temperature(double temperature);

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
    double temperature() const
    {
        return temperature_;
    }

    void bind_textures() const
    {
        brownian_wrapper<dimension, float_type, rng_type>::param.bind(g_param_);
    }

private:
    typedef typename particle_type::position_array_type position_array_type;
    typedef typename particle_type::image_array_type image_array_type;
    typedef typename particle_type::velocity_array_type velocity_array_type;
    typedef typename particle_type::force_array_type force_array_type;
    typedef typename particle_type::torque_array_type torque_array_type;

    typedef utility::profiler::accumulator_type accumulator_type;
    typedef utility::profiler::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type integrate;
    };

    /** system state */
    std::shared_ptr<particle_type> particle_;
    /** random number generator */
    std::shared_ptr<random_type> random_;
    /** simulation domain */
    std::shared_ptr<box_type const> box_;
    /** integration time-step */
    float_type timestep_;
    /** temperature of the heat bath */
    float_type temperature_;
    /** diffusion constant */
    matrix_type diff_const_;
    /** diffusion parameters at CUDA device */
    cuda::vector<float4> g_param_;
    /** module logger */
    std::shared_ptr<logger> logger_;
    /** profiling runtime accumulators */
    runtime runtime_;
};

} // namespace integrators
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_INTEGRATORS_BROWNIAN_HPP */
