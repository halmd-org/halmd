/*
 * Copyright © 2008-2011  Peter Colberg
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

    /** (cutoff lengths + neighbour list skin)² */
    cuda::texture<float> rr_cut_skin;
    /** neighbour list length */
    cuda::symbol<unsigned int> neighbour_size;
    /** neighbour list stride */
    cuda::symbol<unsigned int> neighbour_stride;
    /** number of particles in simulation box */
    cuda::symbol<unsigned int> nbox;
    /** positions, tags */
    cuda::texture<float4> r;
    /** update neighbour lists */
    cuda::function<void (
        int*, unsigned int*, unsigned int const*
      , unsigned int, unsigned int
      , cell_size_type
      , vector_type
    )> update_neighbours;

    static from_binning_wrapper kernel;
};

template <int dimension>
from_binning_wrapper<dimension> const& get_from_binning_kernel()
{
    return from_binning_wrapper<dimension>::kernel;
}

} // namespace neighbours
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_NEIGHBOURS_FROM_BINNING_KERNEL_HPP */
