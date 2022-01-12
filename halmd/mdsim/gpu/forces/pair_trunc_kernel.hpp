/*
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
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

#ifndef HALMD_MDSIM_GPU_FORCES_PAIR_TRUNC_KERNEL_HPP
#define HALMD_MDSIM_GPU_FORCES_PAIR_TRUNC_KERNEL_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <halmd/mdsim/type_traits.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace forces {

template <int dimension, typename potential_type>
struct pair_trunc_wrapper
{
    typedef fixed_vector<float, dimension> vector_type;
    typedef typename type_traits<dimension, float>::gpu::coalesced_vector_type coalesced_vector_type;
    /** compute forces using 32 threads per particle and transposed neighbour lists */
    typedef cuda::function<void (
        potential_type
      , float4 const*
      , cudaTextureObject_t // positions, types
      , coalesced_vector_type*
      , unsigned int const*
      , unsigned int
      , float*
      , float*
      , unsigned int
      , unsigned int
      , vector_type
      , bool
      , float
    )> compute_kernel_unroll_force_loop_type;
    /** computer forced with one thread per particle */
    typedef cuda::function<void (
        potential_type
      , float4 const*
      , cudaTextureObject_t // positions, types
      , coalesced_vector_type*
      , unsigned int const*
      , unsigned int
      , unsigned int
      , float*
      , float*
      , unsigned int
      , unsigned int
      , vector_type
      , bool
      , float
    )> compute_kernel_type;

    unsigned int const nparallel_particles;

    /** compute forces only */
    compute_kernel_unroll_force_loop_type compute_unroll_force_loop;
    /** compute forces and auxiliary stuff: internal energy, potential part of stress tensor, ... */
    compute_kernel_unroll_force_loop_type compute_aux_unroll_force_loop;
    /** compute forces only */
    compute_kernel_type compute;
    /** compute forces and auxiliary stuff: internal energy, potential part of stress tensor, ... */
    compute_kernel_type compute_aux;

    static pair_trunc_wrapper kernel;
};

} // namespace mdsim
} // namespace gpu
} // namespace forces
} // namespace halmd

#endif  /* ! HALMD_MDSIM_GPU_FORCES_PAIR_TRUNC_KERNEL_HPP */
