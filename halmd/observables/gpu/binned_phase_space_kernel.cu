/*
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
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

#include <halmd/observables/gpu/binned_phase_space_kernel.hpp>
#include <halmd/utility/gpu/thread.cuh>

namespace halmd {
namespace observables {
namespace gpu {
namespace binned_phase_space_kernel {

/**
 * compute cell indices from particle positions
 */
template <typename coalesced_vector_type, typename vector_type, typename cell_size_type>
__global__ void compute_cell_index(
    coalesced_vector_type const* g_r
  , unsigned int* g_cell_index
  , vector_type cell_length
  , vector_type cell_origin
  , cell_size_type ncell
  , unsigned int npart
)
{
    enum { dimension = vector_type::static_size };

    unsigned int const i = GTID;
    // need to branch here as g_r has no excess particles
    vector_type r = (i < npart) ? static_cast<vector_type>(g_r[i]) : 0;
    r -= cell_origin;

    // compute cell index from position
    //
    // Since the positions are extended beyond the periodic simulation box,
    // we have to account for modulo operations on large negative values;
    // in C99, "a % n" yields values in (-n, n)
    typedef fixed_vector<int, dimension> cell_diff_type;
    cell_diff_type index = element_mod(
        static_cast<cell_diff_type>(floor(element_div(r, cell_length)))
      , static_cast<cell_diff_type>(ncell)
    );
    for (int j = 0; j < dimension; ++j) {
        if (index[j] < 0) { index[j] += ncell[j]; }
    }

    // convert d-dimensional index to 1-dimensional offset
    // use storage order of C-arrays (last dimension is stored contiguously)
    // FIXME check PTX to ensure CUDA unrolls this loop
    unsigned int offset = index[0];
    for (int j = 1; j < dimension; ++j) {
        offset *= ncell[j];
        offset += index[j];
    }
    g_cell_index[i] = offset;
}

/**
 * compute global cell offsets in particle list
 */
__global__ void find_cell_offset(unsigned int* g_cell_index, unsigned int* g_cell_offset, unsigned int npart)
{
    const unsigned int j = g_cell_index[GTID];
    const unsigned int k = (GTID > 0 && GTID < npart) ? g_cell_index[GTID - 1] : j;

    if (GTID == 0 || k < j) {
        // particle marks the start of a cell
        g_cell_offset[j] = GTID;
    }
}

/**
 * assign particles to cells
 */
__global__ void assign_cells(
    int* g_ret
  , unsigned int const* g_cell_index
  , unsigned int const* g_cell_offset
  , unsigned int const* g_itag
  , unsigned int* g_otag
  , unsigned int npart
  , unsigned int ncell
)
{
    __shared__ unsigned int s_offset[1];

    const unsigned int cell = BID;

    if (TID == 0) {
        s_offset[0] = g_cell_offset[cell];
    }
    __syncthreads();

    // global offset of first particle in this block's cell
    const unsigned int offset = s_offset[0];

    for (unsigned int i = TID; i < ncell; i += TDIM) {
        // global offset of this thread's particle
        const unsigned int n = offset + i;
        // mark as virtual particle
        unsigned int tag = PLACEHOLDER;
        // mark as real particle if appropriate
        if (offset != PLACEHOLDER && n < npart && g_cell_index[n] == cell) {
            tag = g_itag[n];
            // return failure if any cell list is fully occupied
            if (i + 1 == ncell) {
                *g_ret = EXIT_FAILURE;
            }
        }
        // store particle index in this block's cell
        g_otag[cell * ncell + i] = tag;
    }
}

/**
 * generate ascending index sequence
 */
__global__ void generate_sequence(unsigned int* g_index, unsigned int npart)
{
    g_index[GTID] = (GTID < npart) ? GTID : -1U;
}

} // namespace binned_phase_space_kernel

template <int dimension>
binned_phase_space_wrapper<dimension> binned_phase_space_wrapper<dimension>::kernel = {
    binned_phase_space_kernel::assign_cells
  , binned_phase_space_kernel::find_cell_offset
  , binned_phase_space_kernel::generate_sequence
  , binned_phase_space_kernel::compute_cell_index
};

template class binned_phase_space_wrapper<3>;
template class binned_phase_space_wrapper<2>;

} // namespace gpu
} // namespace observables
} // namespace halmd
