/*
 * Copyright © 2017 Daniel Kirchner
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

#include <cuda_wrapper/cuda_wrapper.hpp>

namespace halmd {

template<typename kernel_type>
void configure_kernel(kernel_type const& k, cuda::config const& default_dim, size_t smem_factor = 0) {
    int blockSize = k.max_block_size();
    if (!blockSize) {
        cuda::configure(default_dim.grid, default_dim.block, smem_factor * default_dim.threads_per_block());
    } else {
        int gridSize = (default_dim.threads() + blockSize - 1) / blockSize;
        if (size_t(gridSize) * size_t(blockSize) == default_dim.threads()) {
            cuda::configure(gridSize, blockSize, smem_factor * blockSize);
        } else {
            // if exact block size does not match choose previous power of two
            blockSize |= (blockSize >> 1);
            blockSize |= (blockSize >> 2);
            blockSize |= (blockSize >> 4);
            blockSize |= (blockSize >> 8);
            blockSize |= (blockSize >> 16);
            blockSize -= (blockSize >> 1);
            int gridSize = (default_dim.threads() + blockSize - 1) / blockSize;
            if (size_t(gridSize) * size_t(blockSize) == default_dim.threads()) {
                cuda::configure(gridSize, blockSize, smem_factor * blockSize);
            } else {
                // if this does not match either choose default dimensions
                cuda::configure(default_dim.grid, default_dim.block, smem_factor * default_dim.threads_per_block());
            }
        }
    }
}

} // namespace halmd

#endif /* ! HALMD_UTILITY_GPU_CONFIGURE_KERNEL_HPP */
