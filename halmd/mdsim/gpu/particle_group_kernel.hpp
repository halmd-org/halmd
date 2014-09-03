/*
 * Copyright © 2014 Nicolas Höft
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
#ifndef HALMD_MDSIM_GPU_PARTICLE_GROUP_KERNEL_HPP
#define HALMD_MDSIM_GPU_PARTICLE_GROUP_KERNEL_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>

#include <halmd/mdsim/type_traits.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {

template <int dimension>
struct particle_group_wrapper
{
    typedef typename type_traits<dimension, float>::gpu::vector_type vector_type;
    typedef typename type_traits<dimension, float>::gpu::coalesced_vector_type aligned_vector_type;

    /** positions, types */
    cuda::texture<float4> r;
    /** minimum image vectors */
    cuda::texture<aligned_vector_type> image;
    /** velocities, masses */
    cuda::texture<float4> v;
    /** copy particle group data (position, species, velocities, masses) into particle instance */
    cuda::function<void (unsigned int const*, float4*, aligned_vector_type*, float4*, unsigned int)> particle_group_to_particle;

    static particle_group_wrapper const kernel;
};

} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_PARTICLE_GROUP_KERNEL_HPP */
