/*
 * Copyright © 2011  Felix Höfling
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

#include <halmd/observables/gpu/fields/density_kernel.hpp>
#include <halmd/observables/gpu/binned_phase_space_kernel.hpp>
#include <halmd/utility/gpu/thread.cuh>

namespace halmd {
namespace observables {
namespace gpu {
namespace fields {
namespace density_kernel {

using binned_phase_space_kernel::PLACEHOLDER;

/**
 * determine particle density per binning cell
 *
 * Use one thread block per binning cell. The particle counter is incremented
 * in shared memory by atomic operations. The result is normalised by the cell
 * volume.
 *
 * Each thread of a warp has its own accumulator in shared memory.
 */
__global__ void compute(
    unsigned int const* g_cell
  , double* g_density // store as double, which corresponds to the layout of fields::density::result_type
  , float cell_volume
  , unsigned int cell_size
)
{
    __shared__ unsigned int s_count[32]; // warpSize

    // initialise counters
    if (TID < warpSize) {
        s_count[TID] = 0;
    }
    __syncthreads();

    // cell index
    unsigned int const cell = BID;

    // loop over cell if #threads is less than cell_size
    for (unsigned int i = TID; i < cell_size; i += TDIM ) {
        // load particle index from cell list
        unsigned int const idx = g_cell[cell * cell_size + i];

        // increment counter in shared memory
        if (idx != PLACEHOLDER) {
            atomicInc(s_count + (i % warpSize), -1U);
        }
    }
    __syncthreads();

    // first half-warp reduces the accumulators
    if (TID < warpSize / 2) {
        // no __syncthreads() necessary after each of the
        // following lines as long as we access the data via
        // a pointer declared as volatile because the 32 threads
        // in each warp execute in lock-step with each other
        volatile unsigned int* smem = s_count;
        smem[TID] += smem[TID + 16]; // FIXME replace by recursive template function for arbitrary warpSize
        smem[TID] += smem[TID + 8];
        smem[TID] += smem[TID + 4];
        smem[TID] += smem[TID + 2];
        smem[TID] += smem[TID + 1];

        // thread #0 writes the result to global memory
        if (TID == 0) {
            g_density[cell] = smem[0] / cell_volume;
        }
    }
}

} // namespace density_kernel

density_wrapper density_wrapper::kernel = {
    density_kernel::compute
};

} // namespace fields
} // namespace gpu
} // namespace observables
} // namespace halmd
