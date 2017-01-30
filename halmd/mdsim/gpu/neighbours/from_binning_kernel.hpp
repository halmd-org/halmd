/*
 * Copyright © 2008-2011 Peter Colberg
 * Copyright © 2014      Nicolas Höft
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
#include <halmd/utility/gpu/texture.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace neighbours {

template <int dimension>
struct from_binning_wrapper
{
    typedef fixed_vector<float, dimension> vector_type;
    typedef fixed_vector<unsigned int, dimension> cell_size_type;

    /** (cutoff lengths + neighbour list skin)² */
    cuda::texture<float> rr_cut_skin;
    /** positions, ids of particle1 */
    cuda::halmd::texture<float4> r1;
    /** positions, ids of particle2 */
    cuda::halmd::texture<float4> r2;

    /** update neighbour lists */
    cuda::function<void (
        int*
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
    )> update_neighbours;

    /** update neighbour lists that uses a 'naive' implementation */
    cuda::function<void (
        int*
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
    )> update_neighbours_naive;

    static from_binning_wrapper kernel;
};

} // namespace neighbours
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_NEIGHBOURS_FROM_BINNING_KERNEL_HPP */
