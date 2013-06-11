/*
 * Copyright © 2011-2013 Felix Höfling
 * Copyright © 2011-2012 Peter Colberg
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

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <lua.hpp>
#include <memory>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/mdsim/gpu/particle_group.hpp>
#include <halmd/mdsim/type_traits.hpp>
#include <halmd/observables/gpu/density_mode_kernel.hpp>
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
    typedef mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef mdsim::gpu::particle_group particle_group_type;
    typedef observables::utility::wavevector<dimension> wavevector_type;
    typedef observables::samples::density_mode<dimension> sample_type;
    typedef logger logger_type;

    density_mode(
        std::shared_ptr<particle_type const> particle
      , std::shared_ptr<particle_group_type> particle_group
      , std::shared_ptr<wavevector_type const> wavevector
      , std::shared_ptr<logger_type> logger = std::make_shared<logger_type>("density_mode")
    );

    /**
     * Compute density modes from particle group.
     */
    std::shared_ptr<sample_type const> acquire();

    /**
     * Returns wavevector instance.
     */
    std::shared_ptr<wavevector_type const> wavevector() const
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
    typedef typename sample_type::mode_array_type mode_array_type;
    typedef typename mode_array_type::value_type mode_type;
    typedef density_mode_wrapper<dimension> wrapper_type;

    /** system state */
    std::shared_ptr<particle_type const> particle_;
    /** particle forces */
    std::shared_ptr<particle_group_type> particle_group_;
    /** wavevector grid */
    std::shared_ptr<wavevector_type const> wavevector_;
    /** logger instance */
    std::shared_ptr<logger_type> logger_;

    /** cached sample with density modes */
    std::shared_ptr<sample_type> rho_sample_;
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
