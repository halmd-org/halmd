/*
 * Copyright © 2013       Nicolas Höft
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

#ifndef HALMD_MDSIM_GPU_FORCES_TABULATED_EXTERNAL_KERNEL_HPP
#define HALMD_MDSIM_GPU_FORCES_TABULATED_EXTERNAL_KERNEL_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <halmd/mdsim/type_traits.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace forces {

template <int dimension, typename interpolation_type>
struct tabulated_external_wrapper
{
    typedef fixed_vector<float, dimension> vector_type;
    typedef typename type_traits<dimension, float>::gpu::coalesced_vector_type coalesced_vector_type;

    /** compute forces only */
    cuda::function<void (
        float4 const*
      , coalesced_vector_type*
      , float*
      , float*
      , unsigned int
      , vector_type
      , float const*
      , interpolation_type const
      , bool
    )> compute;

    /** compute forces and auxiliary stuff: internal energy, potential part of stress tensor, ... */
    cuda::function<void (
        float4 const*
      , coalesced_vector_type*
      , float*
      , float*
      , unsigned int
      , vector_type
      , float const*
      , interpolation_type const
      , bool
    )> compute_aux;

    static tabulated_external_wrapper const kernel;
};

} // namespace mdsim
} // namespace gpu
} // namespace forces
} // namespace halmd

#endif  /* ! HALMD_MDSIM_GPU_FORCES_TABULATED_EXTERNAL_KERNEL_HPP */
