/*
 * Copyright © 2011  Michael Kopp
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

#ifndef HALMD_MDSIM_GPU_INTEGRATORS_EULER_HPP
#define HALMD_MDSIM_GPU_INTEGRATORS_EULER_HPP

#include <lua.hpp>
#include <memory>

#include <cuda_wrapper/cuda_wrapper.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/gpu/integrators/euler_kernel.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace integrators {

template <int dimension, typename float_type>
class euler
{
public:
    typedef gpu::particle<dimension, float_type> particle_type;
    typedef mdsim::box<dimension> box_type;
    typedef typename particle_type::vector_type vector_type;
    typedef euler_wrapper<float_type, dimension> wrapper_type;

    static void luaopen(lua_State* L);

    euler(
        std::shared_ptr<particle_type> particle
      , std::shared_ptr<box_type const> box
      , double timestep
      , std::shared_ptr<halmd::logger> logger = std::make_shared<halmd::logger>()
    );

    void integrate();
    void set_timestep(double timestep);

    //! returns integration time step
    double timestep() const
    {
        return timestep_;
    }

private:
    typedef typename particle_type::position_array_type position_array_type;
    typedef typename particle_type::image_array_type image_array_type;
    typedef typename particle_type::velocity_array_type velocity_array_type;

    typedef utility::profiler::accumulator_type accumulator_type;
    typedef utility::profiler::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type integrate;
    };

    std::shared_ptr<particle_type> particle_;
    std::shared_ptr<box_type const> box_;
    /** module logger */
    std::shared_ptr<logger> logger_;
    /** integration time-step */
    float_type timestep_;
    /** profiling runtime accumulators */
    runtime runtime_;
};

} // namespace integrators
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_INTEGRATORS_EULER_HPP */
