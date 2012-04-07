/*
 * Copyright © 2011-2012  Felix Höfling and Peter Colberg
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

#ifndef HALMD_OBSERVABLES_GPU_DENSITY_MODE_HPP
#define HALMD_OBSERVABLES_GPU_DENSITY_MODE_HPP

#include <boost/make_shared.hpp>
#include <cuda_wrapper/cuda_wrapper.hpp>
#include <lua.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/clock.hpp>
#include <halmd/mdsim/type_traits.hpp>
#include <halmd/observables/gpu/density_mode_kernel.hpp>
#include <halmd/observables/gpu/samples/phase_space.hpp>
#include <halmd/observables/samples/density_mode.hpp>
#include <halmd/observables/utility/wavevector.hpp>
#include <halmd/utility/profiler.hpp>

namespace halmd {
namespace observables {
namespace gpu {

/**
 * Compute Fourier modes of the particle density.
 *
 * @f$ \rho_{\vec q} = \sum_{i=1}^N \exp(\textrm{i}\vec q \cdot \vec r_i) @f$
 */
template <int dimension, typename float_type>
class density_mode
{
public:
    typedef gpu::samples::phase_space<dimension, float_type> phase_space_type;
    typedef observables::utility::wavevector<dimension> wavevector_type;
    typedef observables::samples::density_mode<dimension> sample_type;
    typedef mdsim::clock clock_type;
    typedef logger logger_type;

    density_mode(
        boost::shared_ptr<wavevector_type const> wavevector
      , boost::shared_ptr<clock_type const> clock
      , boost::shared_ptr<logger_type> logger = boost::make_shared<logger_type>()
    );

    /**
     * Compute density modes from phase space sample and store with given time stamp (simulation step).
     *
     * FIXME operate on unsorted particle_group instead of phase_space
     */
    boost::shared_ptr<sample_type const> acquire(phase_space_type const& phase_space);

    /**
     * Returns wavevector instance.
     */
    boost::shared_ptr<wavevector_type const> wavevector() const
    {
        return wavevector_;
    }

    /**
     * Bind class to Lua.
     */
    static void luaopen(lua_State* L);

private:
    typedef typename mdsim::type_traits<dimension, float>::vector_type vector_type;
    typedef typename mdsim::type_traits<dimension, float>::gpu::coalesced_vector_type gpu_vector_type;
    typedef typename sample_type::mode_vector_type mode_vector_type;
    typedef typename mode_vector_type::value_type mode_type;
    typedef density_mode_wrapper<dimension> wrapper_type;

    /** cached sample with density modes */
    boost::shared_ptr<sample_type> rho_sample_;
    /** wavevector grid */
    boost::shared_ptr<wavevector_type const> wavevector_;
    /** simulation clock */
    boost::shared_ptr<clock_type const> clock_;
    /** logger instance */
    boost::shared_ptr<logger_type> logger_;

    /** total number of wavevectors */
    unsigned int nq_;
    /** grid and block dimensions for CUDA calls */
    cuda::config const dim_;
    /** wavevectors */
    cuda::vector<gpu_vector_type> g_q_;
    /** block sums of sin(q r) for each wavevector on the device */
    cuda::vector<float> g_sin_block_;
    /** block sums of cos(q r) for each wavevector on the device */
    cuda::vector<float> g_cos_block_;
    /** sin(q r) for each wavevector on the device */
    cuda::vector<float> g_sin_;
    /** cos(q r) for each wavevector on the device */
    cuda::vector<float> g_cos_;
    /** sin(q r) for each wavevector as page-locked host memory */
    cuda::host::vector<float> h_sin_;
    /** cos(q r) for each wavevector as page-locked host memory */
    cuda::host::vector<float> h_cos_;

    typedef halmd::utility::profiler profiler_type;
    typedef typename profiler_type::accumulator_type accumulator_type;
    typedef typename profiler_type::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type acquire;
    };

    /** profiling runtime accumulators */
    runtime runtime_;
};

} // namespace observables
} // namespace gpu
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_GPU_DENSITY_MODE_HPP */
