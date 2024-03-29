/*
 * Copyright © 2021      Jaslo Ziska
 * Copyright © 2014      Nicolas Höft
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

#ifndef HALMD_MDSIM_GPU_NEIGHBOURS_FROM_BINNING_KERNEL_HPP
#define HALMD_MDSIM_GPU_NEIGHBOURS_FROM_BINNING_KERNEL_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace neighbours {

template <int dimension>
struct from_binning_wrapper
{
    typedef fixed_vector<float, dimension> vector_type;
    typedef fixed_vector<unsigned int, dimension> cell_size_type;

    /** update neighbour lists */
    typedef cuda::function<void (
        cudaTextureObject_t // (cutoff distances + neighbour list skin)²
      , cudaTextureObject_t // positions, IDs of particle1
      , cudaTextureObject_t // positions, IDs of particle2
      , int*
      , unsigned int*
      , unsigned int
      , unsigned int
      , unsigned int const*
      , unsigned int const*
      , unsigned int
      , unsigned int
      , unsigned int
      , cell_size_type
      , vector_type
    )> update_neighbours_function_type;

    /** update neighbour lists that uses a 'naive' implementation */
    typedef cuda::function<void (
        cudaTextureObject_t // (cutoff distances + neighbour list skin)²
      , cudaTextureObject_t // positions, IDs of particle2
      , int*
      , float4 const*
      , unsigned int
      , bool
      , unsigned int*
      , unsigned int
      , unsigned int
      , unsigned int const*
      , unsigned int
      , unsigned int
      , cell_size_type
      , vector_type
      , unsigned int
      , vector_type
    )> update_neighbours_naive_function_type;

    struct functions
    {
        update_neighbours_function_type update_neighbours;
        update_neighbours_naive_function_type update_neighbours_naive;
    };

    functions unroll_force_loop;
    functions normal;

    static from_binning_wrapper kernel;
};

} // namespace neighbours
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_NEIGHBOURS_FROM_BINNING_KERNEL_HPP */
