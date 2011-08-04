/*
 * Copyright © 2008-2010  Peter Colberg
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

#ifndef HALMD_MDSIM_GPU_SORTS_HILBERT_KERNEL_HPP
#define HALMD_MDSIM_GPU_SORTS_HILBERT_KERNEL_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>

#include <halmd/mdsim/type_traits.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace sorts {

template <int dimension>
struct hilbert_wrapper
{
    typedef typename type_traits<dimension, float>::gpu::vector_type builtin_vector_type;
    typedef typename type_traits<dimension, float>::gpu::coalesced_vector_type aligned_vector_type;

    /** Hilbert space-filling curve recursion depth */
    cuda::symbol<unsigned int> depth;
    /** cubic box edgle length */
    cuda::symbol<builtin_vector_type> box_length;
    /** positions, types */
    cuda::texture<float4> r;
    /** minimum image vectors */
    cuda::texture<aligned_vector_type> image;
    /** velocities, tags */
    cuda::texture<float4> v;
    /** generate Hilbert space-filling curve */
    cuda::function<void (float4 const*, unsigned int*)> map;
    /** generate ascending index sequence */
    cuda::function<void (unsigned int*)> gen_index;
    /** order particles after given permutation */
    cuda::function<void (unsigned int const*, float4*, aligned_vector_type*, float4*, unsigned int*)> order_particles;

    static hilbert_wrapper const kernel;
};

template <int dimension>
hilbert_wrapper<dimension> const& get_hilbert_kernel()
{
    return hilbert_wrapper<dimension>::kernel;
}

} // namespace mdsim
} // namespace gpu
} // namespace sorts
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_SORTS_HILBERT_KERNEL_HPP */
