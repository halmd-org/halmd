/*
 * Copyright © 2008-2011  Felix Höfling and Peter Colberg
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

#ifndef HALMD_OBSERVABLES_GPU_BINNED_PHASE_SPACE_KERNEL_HPP
#define HALMD_OBSERVABLES_GPU_BINNED_PHASE_SPACE_KERNEL_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <halmd/mdsim/type_traits.hpp>

namespace halmd {
namespace observables {
namespace gpu {

namespace binned_phase_space_kernel {

/** placeholder particle */
enum { PLACEHOLDER = -1U };

} // namespace binned_phase_space_kernel

template <int dimension>
struct binned_phase_space_wrapper
{
    typedef typename mdsim::type_traits<dimension, float>::vector_type vector_type;
    typedef typename mdsim::type_traits<dimension, float>::gpu::coalesced_vector_type coalesced_vector_type;
    typedef typename mdsim::type_traits<dimension, unsigned int>::vector_type cell_index_type;

    /** assign particles to cells */
    cuda::function<void (
        int*
      , unsigned int const*
      , unsigned int const*
      , unsigned int const*
      , unsigned int*
      , unsigned int
      , unsigned int
    )> assign_cells;
    /** compute global cell offsets in particle list */
    cuda::function<void (unsigned int*, unsigned int*, unsigned int)> find_cell_offset;
    /** generate ascending index sequence */
    cuda::function<void (unsigned int*, unsigned int)> generate_sequence;
    /** compute cell indices for particle positions */
    cuda::function<void (
        coalesced_vector_type const*
      , unsigned int*
      , vector_type
      , vector_type
      , cell_index_type
      , unsigned int
    )> compute_cell_index;

    static binned_phase_space_wrapper kernel;
};

} // namespace gpu
} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_GPU_BINNED_PHASE_SPACE_KERNEL_HPP */
