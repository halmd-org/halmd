/*
 * Copyright © 2008-2010  Peter Colberg
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

#ifndef HALMD_OBSERVABLES_GPU_KINETIC_ENERGY_DENSITY_MODE_KERNEL_CUH
#define HALMD_OBSERVABLES_GPU_KINETIC_ENERGY_DENSITY_MODE_KERNEL_CUH

#include <cuda_wrapper/cuda_wrapper.hpp>

#include <halmd/mdsim/type_traits.hpp>

namespace halmd {
namespace observables {
namespace gpu {

template <int dimension>
struct kinetic_energy_density_mode_wrapper
{
    typedef typename mdsim::type_traits<dimension, float>::gpu::vector_type vector_type;
    typedef typename mdsim::type_traits<dimension, float>::gpu::coalesced_vector_type coalesced_vector_type;

    /** compute density modes of kinetic energy for all particles of a single species */
    cuda::function<void (
        cudaTextureObject_t
      , float4 const*
      , float4 const*
      , unsigned int const*, int
      , float*, float*, int
    )> compute;
    /** finalise computation by summing block sums per wavevector */
    cuda::function<void (float const*, float const*, float*, float*, int, int)> finalise;

    static kinetic_energy_density_mode_wrapper kernel;
};

template <int dimension>
kinetic_energy_density_mode_wrapper<dimension>& get_kinetic_energy_density_mode_kernel()
{
    return kinetic_energy_density_mode_wrapper<dimension>::kernel;
}

} // namespace observables
} // namespace gpu
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_GPU_KINETIC_ENERGY_DENSITY_MODE_KERNEL_CUH */
