/*
 * Copyright © 2014      Felix Höfling
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

#include <halmd/mdsim/gpu/binning_kernel.hpp>
#include <halmd/mdsim/gpu/particle_kernel.cuh>
#include <halmd/utility/gpu/thread.cuh>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace binning_kernel {

/**
 * compute cell indices for given particle positions
 */
template <typename vector_type, typename cell_size_type>
inline __device__ unsigned int compute_cell_index(
    vector_type r
  , vector_type cell_length
  , cell_size_type ncell
)
{
    enum { dimension = vector_type::static_size };

    cell_size_type index = element_mod(
        static_cast<cell_size_type>(element_div(r, cell_length) + static_cast<vector_type>(ncell))
      , ncell
    );
    // FIXME check PTX to ensure CUDA unrolls this loop
    unsigned int offset = index[dimension - 1];
    for (int i = dimension - 2; i >= 0; i--) {
        offset *= ncell[i];
        offset += index[i];
    }
    return offset;
}

/**
 * compute cell indices for particle positions
 */
template <unsigned int dimension>
__global__ void compute_cell(
    float4 const* g_r
  , unsigned int* g_cell
  , fixed_vector<float, dimension> cell_length
  , fixed_vector<unsigned int, dimension> ncell
)
{
    fixed_vector<float, dimension> r;
    unsigned int type;
    tie(r, type) <<= g_r[GTID];
    if (type == -1U) {
        g_cell[GTID] = -1U;
    } else {
        g_cell[GTID] = compute_cell_index(r, cell_length, ncell);
    }
}

/**
 * compute global cell offsets in particle list
 */
__global__ void find_cell_offset(unsigned int* g_cell, unsigned int* g_cell_offset, unsigned int const nbox)
{
    const unsigned int j = g_cell[GTID];
    const unsigned int k = (GTID > 0 && GTID < nbox) ? g_cell[GTID - 1] : j;

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
  , unsigned int const* g_cell
  , unsigned int const* g_cell_offset
  , unsigned int const* g_iid
  , unsigned int* g_oid
  , unsigned int const nbox
  , unsigned int const cell_size
)
{
    __shared__ unsigned int s_offset[1];

    if (TID == 0) {
        s_offset[0] = g_cell_offset[BID];
    }
    __syncthreads();

    // global offset of first particle in this block's cell
    const unsigned int offset = s_offset[0];

    // iterate over cell list
    for (unsigned int i = TID; i < cell_size; i += TDIM) {
        // global offset of this thread's particle
        const unsigned int n = offset + i;
        // mark as virtual particle
        unsigned int id = particle_kernel::placeholder;
        // mark as real particle if appropriate
        if (offset != particle_kernel::placeholder && n < nbox && g_cell[n] == BID) {
            id = g_iid[n];
        }
        // return failure if any cell list is fully occupied
        if (id != particle_kernel::placeholder && (i + 1) == cell_size) {
            *g_ret = EXIT_FAILURE;
        }
        // store particle in this block's cell
        g_oid[BID * cell_size + i] = id;
    }
}

/**
 * generate ascending index sequence
 */
__global__ void gen_index(unsigned int* g_index, unsigned int const nbox)
{
    g_index[GTID] = (GTID < nbox) ? GTID : 0;
}

} // namespace binning_kernel

template <int dimension>
binning_wrapper<dimension> binning_wrapper<dimension>::kernel = {
    binning_kernel::assign_cells
  , binning_kernel::find_cell_offset
  , binning_kernel::gen_index
  , binning_kernel::compute_cell<dimension>
};

template class binning_wrapper<3>;
template class binning_wrapper<2>;

} // namespace gpu
} // namespace mdsim
} // namespace halmd
