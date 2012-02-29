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

#ifndef HALMD_OBSERVABLES_GPU_PHASE_SPACE_HPP
#define HALMD_OBSERVABLES_GPU_PHASE_SPACE_HPP

#include <lua.hpp>

#include <halmd/mdsim/box.hpp>
#include <halmd/mdsim/clock.hpp>
#include <halmd/observables/gpu/samples/particle_group.hpp>
#include <halmd/observables/gpu/samples/phase_space.hpp>
#include <halmd/observables/host/samples/phase_space.hpp>
#include <halmd/observables/phase_space.hpp>
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
  : public observables::phase_space<dimension>
{
public:
    typedef observables::phase_space<dimension> _Base;
    typedef gpu::samples::phase_space<dimension, float_type> sample_type;
    typedef gpu::samples::particle_group<dimension, float_type> particle_group_type;
    typedef typename particle_group_type::particle_type particle_type;
    typedef mdsim::box<dimension> box_type;
    typedef mdsim::clock clock_type;
    typedef logger logger_type;
    typedef typename sample_type::vector_type vector_type;

    static void luaopen(lua_State* L);

    // construct from particle group
    phase_space(
        boost::shared_ptr<sample_type> sample
      , boost::shared_ptr<particle_group_type const> particle_group
      , boost::shared_ptr<box_type const> box
      , boost::shared_ptr<clock_type const> clock
      , boost::shared_ptr<logger_type> logger = boost::make_shared<logger_type>()
    );

    virtual void acquire();

private:
    typedef halmd::utility::profiler profiler_type;
    typedef typename profiler_type::accumulator_type accumulator_type;
    typedef typename profiler_type::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type acquire;
        accumulator_type reset;
    };

private:
    boost::shared_ptr<sample_type> sample_;
    boost::shared_ptr<particle_group_type const> particle_group_;
    boost::shared_ptr<box_type const> box_;
    boost::shared_ptr<clock_type const> clock_;
    boost::shared_ptr<logger_type> logger_;
    /** profiling runtime accumulators */
    runtime runtime_;
};

/**
 * Sample phase_space from GPU memory to host memory
 */
template <int dimension, typename float_type>
class phase_space<host::samples::phase_space<dimension, float_type> >
  : public observables::phase_space<dimension>
{
public:
    typedef observables::phase_space<dimension> _Base;
    typedef host::samples::phase_space<dimension, float_type> sample_type;
    typedef gpu::samples::particle_group<dimension, float_type> particle_group_type;
    typedef typename particle_group_type::particle_type particle_type;
    typedef mdsim::box<dimension> box_type;
    typedef fixed_vector<float_type, dimension> vector_type;
    typedef logger logger_type;
    typedef mdsim::clock clock_type;

    static void luaopen(lua_State* L);

    // construct from particle group
    phase_space(
        boost::shared_ptr<sample_type> sample
      , boost::shared_ptr<particle_group_type> particle_group
      , boost::shared_ptr<box_type const> box
      , boost::shared_ptr<clock_type const> clock
      , boost::shared_ptr<logger_type> logger = boost::make_shared<logger_type>()
    );

    virtual void acquire();

private:
    typedef halmd::utility::profiler profiler_type;
    typedef typename profiler_type::accumulator_type accumulator_type;
    typedef typename profiler_type::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type acquire;
        accumulator_type reset;
    };

    boost::shared_ptr<sample_type> sample_;
    boost::shared_ptr<particle_group_type> particle_group_;
    boost::shared_ptr<box_type const> box_;
    boost::shared_ptr<clock_type const> clock_;
    boost::shared_ptr<logger_type> logger_;

    /** buffer data from gpu::particle in page-locked host memory */
    cuda::host::vector<float4> h_r_;
    cuda::host::vector<typename particle_type::gpu_vector_type> h_image_;
    cuda::host::vector<float4> h_v_;

    /** profiling runtime accumulators */
    runtime runtime_;
};

} // namespace observables
} // namespace gpu
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_GPU_PHASE_SPACE_HPP */
