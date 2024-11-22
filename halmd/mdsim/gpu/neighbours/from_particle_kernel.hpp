/*
 * Copyright © 2021      Jaslo Ziska
 * Copyright © 2008-2011 Peter Colberg
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

#ifndef HALMD_MDSIM_GPU_NEIGHBOURS_FROM_PARTICLE_KERNEL_HPP
#define HALMD_MDSIM_GPU_NEIGHBOURS_FROM_PARTICLE_KERNEL_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <halmd/mdsim/type_traits.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace neighbours {

template <int dimension>
struct from_particle_wrapper
{
    /** update neighbour lists */
    typedef cuda::function<void (
        cudaTextureObject_t // (cutoff distances + neighbour list skin)²
      , float4 const*
      , unsigned int
      , float4 const*
      , unsigned int
      , unsigned int
      , unsigned int
      , fixed_vector<float, dimension>
      , unsigned int*
      , unsigned int
      , unsigned int
      , int*
    )> update_function_type;

    update_function_type update_unroll_force_loop;
    update_function_type update;

    static from_particle_wrapper kernel;
};

template <int dimension>
from_particle_wrapper<dimension>& get_from_particle_kernel()
{
    return from_particle_wrapper<dimension>::kernel;
}

} // namespace neighbours
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_NEIGHBOURS_FROM_PARTICLE_KERNEL_HPP */
