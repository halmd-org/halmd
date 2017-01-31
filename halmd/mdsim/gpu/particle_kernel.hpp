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
#include <halmd/utility/gpu/texture.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {

template <typename float_type, int dimension>
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
    cuda::halmd::texture<aligned_hp_vector_type> r;
    /** minimum image vectors */
    cuda::texture<aligned_vector_type> image;
    /** velocities, masses */
    cuda::halmd::texture<aligned_hp_vector_type> v;
    /** IDs */
    cuda::texture<unsigned int> id;
    /** initialize particle positions and species, velocity and mass */
    cuda::function<void (ptr_type, ptr_type, unsigned int)> initialize;
    /** rearrange particles by a given permutation */
    cuda::function<void (unsigned int const*, ptr_type, aligned_vector_type*, ptr_type, unsigned int*, unsigned int)> rearrange;
    static particle_wrapper const kernel;
};

template <typename float_type, int dimension>
particle_wrapper<float_type, dimension> const& get_particle_kernel()
{
    return particle_wrapper<float_type, dimension>::kernel;
}

} // namespace mdsim
} // namespace gpu
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_PARTICLE_KERNEL_HPP */
