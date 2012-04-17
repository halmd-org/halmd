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

#ifndef HALMD_MDSIM_GPU_BINNING_KERNEL_HPP
#define HALMD_MDSIM_GPU_BINNING_KERNEL_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <halmd/mdsim/type_traits.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {

template <int dimension>
struct binning_wrapper
{
    typedef fixed_vector<float, dimension> vector_type;
    typedef fixed_vector<unsigned int, dimension> index_type;

    /** number of particles in simulation box */
    cuda::symbol<unsigned int> nbox;
    /** assign particles to cells */
    cuda::function<void (int*, unsigned int const*, unsigned int const*, unsigned int const*, unsigned int*)> assign_cells;
    /** compute global cell offsets in particle list */
    cuda::function<void (unsigned int*, unsigned int*)> find_cell_offset;
    /** generate ascending index sequence */
    cuda::function<void (unsigned int*)> gen_index;
    /** compute cell indices for particle positions */
    cuda::function<void (float4 const*, unsigned int*, vector_type, index_type)> compute_cell;

    static binning_wrapper kernel;
};

template <int dimension>
binning_wrapper<dimension> const& get_binning_kernel()
{
    return binning_wrapper<dimension>::kernel;
}

} // namespace mdsim
} // namespace gpu
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_BINNING_KERNEL_HPP */
