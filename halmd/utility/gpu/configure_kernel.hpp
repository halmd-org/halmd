/*
 * Copyright © 2020-2021 Jaslo Ziska
 * Copyright © 2017-2021 Felix Höfling
 * Copyright © 2017      Daniel Kirchner
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

#ifndef HALMD_UTILITY_GPU_CONFIGURE_KERNEL_HPP
#define HALMD_UTILITY_GPU_CONFIGURE_KERNEL_HPP

#include <halmd/io/logger.hpp>
#include <halmd/utility/gpu/device.hpp>

#include <algorithm>
#include <cmath>
#include <cuda_wrapper/cuda_wrapper.hpp>

namespace halmd {

/**
 * Find optimal block size for a given CUDA kernel for a fixed or minimal total
 * number of threads.
 *
 * @param k: CUDA kernel, instance of cuda::function from cuda_wrapper
 * @param dim: initial configuration dimensions, specifies the total number of threads
 * @param fixed_total_threads: enforce the total number of threads given by dim
 * @param smem_per_thread: shared memory requirements per thread (in bytes)
 *
 * @returns: configuration parameters
 *
 */
template <typename kernel_type>
cuda::config compute_kernel_dimensions(
    kernel_type& k
  , cuda::config dim // pass by copy so that we can change it
  , bool fixed_total_threads
  , size_t smem_per_thread = 0
)
{
    cuda::device::properties prop(device::get());
    size_t warp_size = prop.warp_size();
    int threads = dim.threads();

    assert(threads % warp_size == 0);
    int warps = threads / warp_size;        // total number of warps in the grid

    int min_grid_size = k.min_grid_size();
    int max_block_size = k.max_block_size();
    if (smem_per_thread > 0) {
        max_block_size = std::min(max_block_size, int(prop.shared_mem_per_block() / smem_per_thread));
        max_block_size = (max_block_size / warp_size) * warp_size; // round down to a multiple of warp_size
    }
    int max_warps_per_block = max_block_size / warp_size;

    if (fixed_total_threads) {
        // If the number of threads is fixed we need to find a number of warps
        // per block which is a divisor of the total number of warps, which
        // does not exceed the maximum number of warps per block and which
        // yields a grid size that is close (but smaller or equal than) the
        // "minimal" grid size or a multiple of it.
        //
        // We use a cost model that assumes that the execution of a single warp
        // does not depend on block size; second, blocks that fit within the
        // minimal grid size run perfectly in parallel. Since this model would
        // yield very small block sizes, we start at the largest possible block
        // size and decrease the size only if the cost is reduced by a
        // significant factor.
        int warps_per_block = dim.threads_per_block() / warp_size;
        int cost = int((dim.blocks_per_grid() + min_grid_size - 1) / min_grid_size) * warps_per_block;
        double factor = 0.8;

        for (int w = max_warps_per_block; w >= 1; --w) {
            if (warps % w == 0 ) {
                // estimate runtime cost
                int grid_size = warps / w;
                int c = int((grid_size + min_grid_size - 1) / min_grid_size) * w;
                // update configuration
                if (c < factor * cost) {
                    cost = c;
                    warps_per_block = w;
                    LOG_TRACE("update: warps: " << w << ", cost: " << c);
                }
            }
        }

        // store obtained execution dimensions
        assert(warps % warps_per_block == 0);
        dim = cuda::config(warps / warps_per_block, warps_per_block * warp_size);
    } else {
        // When the number of threads is not fixed, optimal execution dimensions can be calculated.
        // The goal is to have the grid size equal to a multiple of min_grid_size or (slightly) below.
        int r = max_warps_per_block * min_grid_size;
        int grid_size = int((warps + r - 1) / r) * min_grid_size;
        int block_size = int((warps + grid_size - 1) / grid_size) * warp_size;
        assert(block_size <= max_block_size);
        assert(grid_size * block_size >= threads);
        dim = cuda::config(grid_size, block_size);
    }

    return dim;
}

/**
 * Find optimal block size for a given CUDA kernel using compute_kernel_dimensions()
 * and configure that kernel.
 *
 * @param k: CUDA kernel, instance of cuda::function from cuda_wrapper
 * @param dim: initial configuration dimensions, specifies the total number of threads
 * @param fixed_total_threads: enforce the total number of threads given by dim
 * @param smem_per_thread: shared memory requirements per thread (in bytes)
 *
 * @returns: configuration parameters
 *
 */
template <typename kernel_type>
cuda::config configure_kernel(
    kernel_type& k
  , cuda::config dim
  , bool fixed_total_threads
  , size_t smem_per_thread = 0
)
{
    dim = compute_kernel_dimensions(k, dim, fixed_total_threads, smem_per_thread);

    LOG_TRACE("Configuring CUDA kernel for " << dim.blocks_per_grid() << " blocks of " << dim.threads_per_block() << " threads each");
    k.configure(dim.grid, dim.block, smem_per_thread * dim.threads_per_block());

    return dim;
}

/**
 * Configure a given CUDA kernel for a given the total number of threads.
 *
 * @param k: CUDA kernel, instance of cuda::function from cuda_wrapper
 * @param total_threads: total number of threads
 * @param fixed_total_threads: enforce the total number of threads
 * @param smem_per_thread: shared memory requirements per thread (in bytes)
 *
 * @returns: configuration parameters
 */
template <typename kernel_type>
cuda::config configure_kernel(
    kernel_type& k
  , size_t total_threads
  , bool fixed_total_threads = false
  , size_t smem_per_thread = 0
)
{
    unsigned int block_size = 32 << DEVICE_SCALE; // default value for old CUDA versions: not too large, and not too small
    unsigned int grid_size = (total_threads + block_size - 1) / block_size;

    // ensure that the number of threads is unchanged when it is fixed
    assert(fixed_total_threads || grid_size * block_size == total_threads);

    return configure_kernel(k, cuda::config(grid_size, block_size), fixed_total_threads, smem_per_thread);
}

} // namespace halmd

#endif /* ! HALMD_UTILITY_GPU_CONFIGURE_KERNEL_HPP */
