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

template <int dimension, typename float_type>
struct particle_wrapper
{
    typedef typename type_traits<dimension, float>::gpu::vector_type vector_type;
    typedef typename type_traits<dimension, float>::gpu::coalesced_vector_type aligned_vector_type;
    typedef typename type_traits<4, float_type>::gpu::coalesced_vector_type aligned_hp_vector_type;
    typedef typename type_traits<4, float_type>::gpu::ptr_type ptr_type;

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
    cuda::function<void (unsigned int const*, ptr_type, aligned_vector_type*, ptr_type, unsigned int*, unsigned int)> rearrange;
    static particle_wrapper const kernel;
};

template <int dimension, typename float_type>
particle_wrapper<dimension, float_type> const& get_particle_kernel()
{
    return particle_wrapper<dimension, float_type>::kernel;
}

template<typename T>
struct particle_initialize_wrapper
{
    cuda::function<void (T*, T, T, unsigned int)> initialize;
    static particle_initialize_wrapper const kernel;
};

#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
template<size_t dimension>
struct dsfloat_particle_initialize_wrapper
{
    typedef typename type_traits<dimension, dsfloat>::gpu::ptr_type ptr_type;
    typedef typename type_traits<dimension, float>::gpu::coalesced_vector_type type;
    cuda::function<void (ptr_type, type, type, unsigned int)> initialize;
    static dsfloat_particle_initialize_wrapper const kernel;
};
#endif

} // namespace mdsim
} // namespace gpu
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_PARTICLE_KERNEL_HPP */
