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

#ifndef HALMD_OBSERVABLES_GPU_PHASE_SPACE_HPP
#define HALMD_OBSERVABLES_GPU_PHASE_SPACE_HPP

#include <lua.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/clock.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/mdsim/gpu/particle_group.hpp>
#include <halmd/observables/gpu/samples/phase_space.hpp>
#include <halmd/observables/host/samples/phase_space.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd {
namespace observables {
namespace gpu {

template <typename sample_type>
class phase_space;

/**
 * Sample phase_space from GPU memory to GPU memory
 */
template <int dimension, typename float_type>
class phase_space<gpu::samples::phase_space<dimension, float_type> >
{
public:
    typedef gpu::samples::phase_space<dimension, float_type> sample_type;
    typedef mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef mdsim::gpu::particle_group particle_group_type;
    typedef mdsim::box<dimension> box_type;
    typedef mdsim::clock clock_type;

    /**
     * Construct phase_space sampler from particle group.
     */
    phase_space(
        std::shared_ptr<particle_type> particle
      , std::shared_ptr<particle_group_type> particle_group
      , std::shared_ptr<box_type const> box
      , std::shared_ptr<clock_type const> clock
      , std::shared_ptr<logger> logger = std::make_shared<logger>()
    );

    /**
     * Acquire phase_space sample.
     */
    std::shared_ptr<sample_type const> acquire();

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    typedef typename particle_type::position_array_type position_array_type;
    typedef typename particle_type::image_array_type image_array_type;
    typedef typename particle_type::velocity_array_type velocity_array_type;
    typedef typename particle_group_type::array_type group_array_type;
    typedef typename sample_type::vector_type vector_type;

    /** particle instance to particle group */
    std::shared_ptr<particle_type> particle_;
    /** particle group */
    std::shared_ptr<particle_group_type> particle_group_;
    /** simulation box */
    std::shared_ptr<box_type const> box_;
    /** simulation clock */
    std::shared_ptr<clock_type const> clock_;
    /** logger instance */
    std::shared_ptr<logger> logger_;
    /** cached phase_space sample */
    std::shared_ptr<sample_type> sample_;

    typedef halmd::utility::profiler::accumulator_type accumulator_type;
    typedef halmd::utility::profiler::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type acquire;
        accumulator_type reset;
    };

    /** profiling runtime accumulators */
    runtime runtime_;
};

/**
 * Sample phase_space from GPU memory to host memory
 */
template <int dimension, typename float_type>
class phase_space<host::samples::phase_space<dimension, float_type> >
{
public:
    typedef host::samples::phase_space<dimension, float_type> sample_type;
    typedef mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef mdsim::gpu::particle_group particle_group_type;
    typedef mdsim::box<dimension> box_type;
    typedef fixed_vector<float_type, dimension> vector_type;
    typedef mdsim::clock clock_type;

    /**
     * Construct phase_space sampler from particle group.
     */
    phase_space(
        std::shared_ptr<particle_type> particle
      , std::shared_ptr<particle_group_type> particle_group
      , std::shared_ptr<box_type const> box
      , std::shared_ptr<clock_type const> clock
      , std::shared_ptr<logger> logger = std::make_shared<logger>()
    );

    /**
     * Acquire phase_space sample.
     */
    std::shared_ptr<sample_type const> acquire();

    /**
     * Set particles from phase_space sample.
     */
    void set(std::shared_ptr<sample_type const> sample);

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    typedef typename particle_type::position_array_type position_array_type;
    typedef typename particle_type::image_array_type image_array_type;
    typedef typename particle_type::velocity_array_type velocity_array_type;
    typedef typename particle_group_type::array_type group_array_type;

    /** particle instance to particle group */
    std::shared_ptr<particle_type> particle_;
    /** particle group */
    std::shared_ptr<particle_group_type> particle_group_;
    /** simulation box */
    std::shared_ptr<box_type const> box_;
    /** simulation clock */
    std::shared_ptr<clock_type const> clock_;
    /** logger instance */
    std::shared_ptr<logger> logger_;
    /** cached phase_space sample */
    std::shared_ptr<sample_type> sample_;

    /** buffered positions in page-locked host memory */
    cuda::host::vector<float4> h_r_;
    /** buffered periodic image vectors in page-locked host memory */
    cuda::host::vector<typename particle_type::gpu_vector_type> h_image_;
    /** buffered velocities in page-locked host memory */
    cuda::host::vector<float4> h_v_;
    /** GPU threads per block */
    unsigned int threads_;

    typedef halmd::utility::profiler::accumulator_type accumulator_type;
    typedef halmd::utility::profiler::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type acquire;
        accumulator_type reset;
        accumulator_type set;
    };

    /** profiling runtime accumulators */
    runtime runtime_;
};

} // namespace observables
} // namespace gpu
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_GPU_PHASE_SPACE_HPP */
