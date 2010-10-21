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

#ifndef HALMD_MDSIM_GPU_FORCES_PAIR_SHORT_RANGED_KERNEL_HPP
#define HALMD_MDSIM_GPU_FORCES_PAIR_SHORT_RANGED_KERNEL_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <halmd/mdsim/type_traits.hpp>

namespace halmd
{
namespace mdsim { namespace gpu { namespace forces
{
namespace pair_short_ranged_kernel
{

//
// index of cutoff in potential parameters, is assumed to be zero
//
enum {
    /** square of cutoff length */
    RR_CUT = 0
};

} // namespace pair_short_ranged_kernel

template <int dimension, typename potential_type>
struct pair_short_ranged_wrapper
{
    typedef typename type_traits<dimension, float>::gpu::coalesced_vector_type coalesced_vector_type;
    typedef typename type_traits<dimension, float>::gpu::vector_type vector_type;
    typedef typename type_traits<dimension, float>::gpu::stress_tensor_type stress_tensor_type;

    /** compute forces, internal energy, and potential part of stress tensor */
    cuda::function<void (coalesced_vector_type*, unsigned int*, float*, stress_tensor_type*)> compute;
    /** cubic box edgle length */
    cuda::symbol<vector_type> box_length;
    /** number of placeholders per neighbour list */
    cuda::symbol<unsigned int> neighbour_size;
    /** neighbour list stride */
    cuda::symbol<unsigned int> neighbour_stride;
    /** positions, types */
    cuda::texture<float4> r;

    static pair_short_ranged_wrapper const kernel;
};

}}} // namespace mdsim::gpu::forces

} // namespace halmd

#endif  /* ! HALMD_MDSIM_GPU_FORCES_PAIR_SHORT_RANGED_KERNEL_HPP */
