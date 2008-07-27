/* Parallel exclusive prefix sum
 *
 * Copyright (C) 2008  Peter Colberg
 *
 * This program is free software: you can redistribute it and/or modify
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

#ifndef MDSIM_GPU_SCAN_GLUE_HPP
#define MDSIM_GPU_SCAN_GLUE_HPP

#include <cuda_wrapper.hpp>

namespace mdsim { namespace gpu { namespace scan
{

enum {
    // number of shared memory banks
    SHMEM_BANKS = 16,
};

/**
 * returns array index with offset for bank conflict free shared memory access
 */
__device__ __host__ inline uint boff(uint const& i)
{
    return i + i / SHMEM_BANKS;
}

extern cuda::function<void (uint const*, uint*, uint*, const uint)> block_prefix_sum;
extern cuda::function<void (uint const*, uint const*, uint*, const uint)> add_block_sums;

}}} // namespace mdsim::gpu::scan

#endif /* ! MDSIM_GPU_SCAN_GLUE_HPP */
