/*
 * Copyright © 2010  Peter Colberg
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

#ifndef HALMD_MDSIM_GPU_PARTICLE_KERNEL_HPP
#define HALMD_MDSIM_GPU_PARTICLE_KERNEL_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>

#include <halmd/mdsim/type_traits.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {

template <int dimension>
struct particle_wrapper
{
    typedef typename type_traits<dimension, float>::gpu::vector_type vector_type;
    typedef typename type_traits<dimension, float>::gpu::coalesced_vector_type aligned_vector_type;

    cuda::symbol<unsigned int> nbox;
    cuda::symbol<unsigned int> ntype;
    cuda::texture<unsigned int> ntypes;
    /** positions, types */
    cuda::texture<float4> r;
    /** minimum image vectors */
    cuda::texture<aligned_vector_type> image;
    /** velocities, masses */
    cuda::texture<float4> v;
    /** IDs */
    cuda::texture<unsigned int> id;
    /** rearrange particles by a given permutation */
    cuda::function<void (unsigned int const*, float4*, aligned_vector_type*, float4*, unsigned int*, unsigned int)> rearrange;
    static particle_wrapper const kernel;
};

template <int dimension>
particle_wrapper<dimension> const& get_particle_kernel()
{
    return particle_wrapper<dimension>::kernel;
}

} // namespace mdsim
} // namespace gpu
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_PARTICLE_KERNEL_HPP */
