/*
 * Copyright © 2020-2021 Jaslo Ziska
 * Copyright © 2017      Daniel Kirchner
 * Copyright © 2017      Felix Höfling
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
#include <cuda_wrapper/cuda_wrapper.hpp>

namespace halmd {

/**
 * Find maximum possible block size for a given CUDA kernel for a fixed or
 * minimal total number of threads and configure that kernel.
 *
 * @param k: CUDA kernel, instance of cuda::function from cuda_wrapper
 * @param default_dim: default configuration dimensions, specifies the total number of threads
 * @param fixed_total_threads: enforce the total number of threads given by default_dim
 * @param smem_per_thread: shared memory requirements per thread (in bytes)
 *
 * @returns: configuration parameters
 *
 */
template <typename kernel_type>
cuda::config configure_kernel(kernel_type& k, cuda::config const& default_dim,
    bool fixed_total_threads, size_t smem_per_thread = 0)
{
    cuda::config dim = default_dim;
    cuda::device::properties prop(device::get());
    size_t warp_size = prop.warp_size();
    int threads = dim.threads();

    assert(threads % warp_size == 0);

    // unsigned long to prevent overflow when calculating warps_per_block
    unsigned long int max_grid_size = prop.max_grid_size().x;
    int min_grid_size = k.min_grid_size();

    int max_block_size = k.max_block_size();
    if (smem_per_thread > 0) {
        max_block_size = std::min(max_block_size, int(prop.shared_mem_per_block() / smem_per_thread));
        max_block_size = (max_block_size / warp_size) * warp_size; // round down to a multiple of warp_size
    }
    int max_warps_per_block = max_block_size / warp_size;

    int warps = threads / warp_size;
    // the (minimum) number of warps_per_block is the number of warps divided by the max_grid_size (rounded up)
    int warps_per_block = (warps + max_grid_size - 1) / max_grid_size;

    if (fixed_total_threads) {
        // If the number of threads is fixed we need to find a number of warps_per_block which is a divisor of the total
        // number of warps and which lies in the interval of maximum and minimum number of warps per block
        while (warps % warps_per_block != 0 &&
               warps_per_block <= max_warps_per_block &&
               warps / warps_per_block >= min_grid_size) {
            warps_per_block++;
        }

        // Only change dim if we found a solution (when warps_per_block is a divisor of the total number of warps).
        // Else we exited the loop because the maximum block size or the minimum grid size were exceeded.
        if (warps % warps_per_block == 0) {
            dim = cuda::config(warps / warps_per_block, warps_per_block * warp_size);
        }
    } else {
        // When the number of threads is not fixed we can just calculate the dimensions with warps_per_block.
        // The grid size is the number of warps divided by the warps_per_block (rounded up).
        // The block size is just the warps_per_thread multiplied by the warp_size.
        dim = cuda::config((warps + warps_per_block - 1) / warps_per_block, warps_per_block * warp_size);
    }

    LOG_TRACE("Configure GPU kernel for " << dim.blocks_per_grid() << " blocks of " << dim.threads_per_block() << " threads each");
    k.configure(dim.grid, dim.block, smem_per_thread * dim.threads_per_block());

    return dim;
}

/**
 * Configure a given CUDA kernel for a given minimal total number of threads.
 *
 * @returns: configuration parameters
 */
template <typename kernel_type>
cuda::config configure_kernel(kernel_type& k, size_t total_threads,
    size_t smem_per_thread = 0)
{
    unsigned int block_size = 32 << DEVICE_SCALE;    // default value for old CUDA versions: not too large, and not too small
    unsigned int grid_size = (total_threads + block_size - 1) / block_size;
    return configure_kernel(k, cuda::config(grid_size, block_size), false, smem_per_thread);
}

} // namespace halmd

#endif /* ! HALMD_UTILITY_GPU_CONFIGURE_KERNEL_HPP */
