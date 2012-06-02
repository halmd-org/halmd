/*
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
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

#ifndef HALMD_MDSIM_GPU_INTEGRATORS_VERLET_HPP
#define HALMD_MDSIM_GPU_INTEGRATORS_VERLET_HPP

#include <lua.hpp>
#include <memory>

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/gpu/integrators/verlet_kernel.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace integrators {

template <int dimension, typename float_type>
class verlet
{
public:
    typedef gpu::particle<dimension, float_type> particle_type;
    typedef mdsim::box<dimension> box_type;
    typedef logger logger_type;
    typedef typename particle_type::vector_type vector_type;

    static void luaopen(lua_State* L);

    verlet(
        std::shared_ptr<particle_type> particle
      , std::shared_ptr<box_type const> box
      , double timestep
      , std::shared_ptr<logger_type> logger = std::make_shared<logger_type>()
    );
    void integrate();
    void finalize();
    void set_timestep(double timestep);

    //! returns integration time-step
    double timestep() const
    {
        return timestep_;
    }

private:
    typedef utility::profiler profiler_type;
    typedef typename profiler_type::accumulator_type accumulator_type;
    typedef typename profiler_type::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type integrate;
        accumulator_type finalize;
    };

    std::shared_ptr<particle_type> particle_;
    std::shared_ptr<box_type const> box_;
    /** module logger */
    std::shared_ptr<logger_type> logger_;
    /** CUDA C++ wrapper */
    verlet_wrapper<dimension> const* wrapper_;
    /** integration time-step */
    float_type timestep_;
    /** profiling runtime accumulators */
    runtime runtime_;
};

} // namespace mdsim
} // namespace gpu
} // namespace integrators
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_INTEGRATORS_VERLET_HPP */
