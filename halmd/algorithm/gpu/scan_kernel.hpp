/*
 * Copyright © 2008-2009  Peter Colberg
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

#ifndef HALMD_ALGORITHM_GPU_SCAN_KERNEL_HPP
#define HALMD_ALGORITHM_GPU_SCAN_KERNEL_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>

namespace halmd {
namespace algorithm {
namespace gpu {

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

template <typename T>
struct scan_wrapper
{
    cuda::function<void (T const*, T*, T*, const uint)> grid_prefix_sum;
    cuda::function<void (T const*, T*, T const*, const uint)> add_block_sums;
    cuda::function<void (T const*, T*, const uint)> block_prefix_sum;
    static scan_wrapper const kernel;
};

template <typename T>
scan_wrapper<T> const& get_scan_kernel()
{
    return scan_wrapper<T>::kernel;
}

} // namespace algorithm
} // namespace gpu
} // namespace halmd

#endif /* ! HALMD_ALGORITHM_GPU_SCAN_KERNEL_HPP */
