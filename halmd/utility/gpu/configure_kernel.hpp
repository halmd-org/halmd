/*
 * Copyright © 2017 Daniel Kirchner
 * Copyright © 2017 Felix Höfling
 * Copyright © 2020 Jaslo Ziska
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

    int block_size = k.max_block_size();
    if (smem_per_thread > 0) {
        size_t shared_mem_per_block = 48 << 10;     // 48 kB
        block_size = std::min(block_size, int(shared_mem_per_block / smem_per_thread));
        block_size = (block_size / 32) * 32; // round block_size down to a multiple of 32 (the warpSize)
    }

    if (block_size > 0) {
        int grid_size = (default_dim.threads() + block_size - 1) / block_size;
        // if the resulting grid size is too small to achieve a good device occupancy
        // use the default dimensions (which has a minimal block size)
        if (grid_size >= k.min_grid_size()) {
            if (!fixed_total_threads || (size_t(grid_size) * size_t(block_size) == default_dim.threads())) {
                dim = cuda::config(grid_size, block_size);
            } else {
                // if exact block size does not match choose previous power of two
                block_size |= (block_size >> 1);
                block_size |= (block_size >> 2);
                block_size |= (block_size >> 4);
                block_size |= (block_size >> 8);
                block_size |= (block_size >> 16);
                block_size -= (block_size >> 1);

                int grid_size = (default_dim.threads() + block_size - 1) / block_size;
                if (size_t(grid_size) * size_t(block_size) == default_dim.threads()) {
                    dim = cuda::config(grid_size, block_size);
                } // otherwise, choose default dimensions
            }
        }
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
