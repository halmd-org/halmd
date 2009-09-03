/* Parallel exclusive prefix sum
 *
 * Copyright © 2008-2009  Peter Colberg
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

#ifndef LJGPU_ALGORITHM_GPU_PREFIX_SUM_HPP
#define LJGPU_ALGORITHM_GPU_PREFIX_SUM_HPP

#include <cuda_wrapper.hpp>
#include <ljgpu/rng/gpu/uint48.cuh>

namespace ljgpu { namespace gpu { namespace prefix_sum
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

extern cuda::function<void (uint const*, uint*, uint*, const uint),
                      void (uint48 const*, uint48*, uint48*, const uint)> grid_prefix_sum;
extern cuda::function<void (uint const*, uint*, uint const*, const uint),
                      void (uint48 const*, uint48*, uint48 const*, const uint)> add_block_sums;
extern cuda::function<void (uint48 const*, uint48*, const uint),
                      void (uint const*, uint*, const uint)> block_prefix_sum;

}}} // namespace ljgpu::gpu::radix_sort

#endif /* ! LJGPU_ALGORITHM_GPU_PREFIX_SUM_HPP */
