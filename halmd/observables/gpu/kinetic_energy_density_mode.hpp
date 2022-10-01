/*
 * Copyright © 2011-2022 Felix Höfling
 * Copyright © 2011-2012 Peter Colberg
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

#ifndef HALMD_OBSERVABLES_GPU_KINETIC_ENERGY_DENSITY_MODE_HPP
#define HALMD_OBSERVABLES_GPU_KINETIC_ENERGY_DENSITY_MODE_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <lua.hpp>
#include <memory>

#include <halmd/io/logger.hpp>
#include <halmd/mdsim/gpu/particle.hpp>
#include <halmd/mdsim/gpu/particle_group.hpp>
#include <halmd/mdsim/type_traits.hpp>
#include <halmd/observables/gpu/kinetic_energy_density_mode_kernel.hpp>
#include <halmd/observables/utility/wavevector.hpp>
#include <halmd/utility/cache.hpp>
#include <halmd/utility/owner_equal.hpp>
#include <halmd/utility/profiler.hpp>
#include <halmd/utility/raw_array.hpp>

namespace halmd {
namespace observables {
namespace gpu {

/**
 * Compute Fourier modes of the kinetic energy density.
 *
 * @f$ \rho_{\vec q} = \sum_{i=1}^N (m_i v_i^2 / 2) \exp(\textrm{i}\vec q \cdot \vec r_i) @f$
 *
 * The result is stored and returned within a std::shared_ptr, allowing
 * efficient copying, e.g., in dynamics::blocking_scheme.  Further, the result
 * may be tracked by std::weak_ptr providing a similar functionality as
 * halmd::cache
 */
template <int dimension, typename float_type>
class kinetic_energy_density_mode
{
public:
    typedef mdsim::gpu::particle<dimension, float_type> particle_type;
    typedef mdsim::gpu::particle_group particle_group_type;
    typedef observables::utility::wavevector<dimension> wavevector_type;
    typedef raw_array<fixed_vector<double, 2>> result_type;

    kinetic_energy_density_mode(
        std::shared_ptr<particle_type const> particle
      , std::shared_ptr<particle_group_type> particle_group
      , std::shared_ptr<wavevector_type const> wavevector
      , std::shared_ptr<halmd::logger> logger = std::make_shared<halmd::logger>("kinetic_energy_density_mode")
    );

    /** Compute current density modes from particle group.
     *
     * The result is re-computed only if the particle positions have been
     * modified. In this case, the data array managed by std::shared_ptr is
     * re-allocated before.
     */
    std::shared_ptr<result_type const> acquire();

    /**
     * Functor wrapping acquire() for class instance stored within std::shared_ptr
     */
    static std::function<std::shared_ptr<result_type const> ()>
    acquisitor(std::shared_ptr<kinetic_energy_density_mode> self)
    {
        return [=]() {
           return self->acquire();
        };
    }

    /**
     * Return wavevector instance passed to constructor
     */
    std::shared_ptr<wavevector_type const> wavevector()
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
    typedef kinetic_energy_density_mode_wrapper<dimension> wrapper_type;

    /** system state */
    std::shared_ptr<particle_type const> particle_;
    /** particle forces */
    std::shared_ptr<particle_group_type> particle_group_;
    /** wavevector grid */
    std::shared_ptr<wavevector_type const> wavevector_;
    /** logger instance */
    std::shared_ptr<logger> logger_;

    /** result for the density modes */
    std::shared_ptr<result_type> result_;
    /** cache observer for particle positions */
    cache<> position_cache_;
    /** cache observer for particle velocities */
    cache<> velocity_cache_;
    /** cache observer for particle group */
    cache<> group_cache_;

    /** total number of wavevectors */
    unsigned int nq_;
    /** grid and block dimensions for CUDA calls */
    cuda::config const dim_;
    /** wavevectors */
    cuda::vector<gpu_vector_type> g_wavevector_;
    /** block sums of sin(q r) for each wavevector on the device */
    cuda::vector<float> g_sin_block_;
    /** block sums of cos(q r) for each wavevector on the device */
    cuda::vector<float> g_cos_block_;
    /** sin(q r) for each wavevector on the device */
    cuda::vector<float> g_sin_;
    /** cos(q r) for each wavevector on the device */
    cuda::vector<float> g_cos_;
    /** sin(q r) for each wavevector as page-locked host memory */
    cuda::memory::host::vector<float> h_sin_;
    /** cos(q r) for each wavevector as page-locked host memory */
    cuda::memory::host::vector<float> h_cos_;

    typedef halmd::utility::profiler::accumulator_type accumulator_type;
    typedef halmd::utility::profiler::scoped_timer_type scoped_timer_type;

    struct runtime
    {
        accumulator_type acquire;
    };

    /** profiling runtime accumulators */
    runtime runtime_;
};

} // namespace gpu
} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_GPU_KINETIC_ENERGY_DENSITY_MODE_HPP */
