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

#ifndef HALMD_ALGORITHM_GPU_SCAN_KERNEL_CUH
#define HALMD_ALGORITHM_GPU_SCAN_KERNEL_CUH

#include <halmd/algorithm/gpu/bits.cuh>
#include <halmd/algorithm/gpu/scan_kernel.hpp>

namespace halmd
{
namespace algorithm { namespace gpu
{
namespace scan_kernel
{

/**
 * blockwise parallel exclusive prefix sum
 */
template <typename T>
__device__ T grid_prefix_sum(T const* g_in, T* g_out, uint const count)
{
    //
    // Prefix Sums and Their Applications,
    // Guy E. Blelloch.
    // CMU-CS-90-190, November 1990.
    //
    // http://www.cs.cmu.edu/~scandal/papers/CMU-CS-90-190.html
    //

    //
    // Parallel Prefix Sum (Scan) with CUDA,
    // Mark Harris, April 2007, NVIDIA Corporation
    //

    extern __shared__ char __s_array[];
    T* const s_array = reinterpret_cast<T*>(__s_array); // work around for CUDA 3.0/3.1
    T block_sum = T(); // value-initialized

    uint const tid = threadIdx.x;
    uint const threads = blockDim.x;
    uint const bid = blockIdx.x;

    // read elements from global memory, or pad with zero
    uint const i1 = 2 * bid * threads + tid;
    uint const i2 = (2 * bid + 1) * threads + tid;
    s_array[boff(tid)] = (i1 < count) ? g_in[i1] : T();
    s_array[boff(threads + tid)] = (i2 < count) ? g_in[i2] : T();
    __syncthreads();

    // up-sweep phase from leaves to root of binary tree
    for (uint d = threads, n = 1; d > 0; d >>= 1, n <<= 1) {
        if (tid < d) {
            s_array[boff(n * (2 * tid + 2) - 1)] += s_array[boff(n * (2 * tid + 1) - 1)];
        }
        __syncthreads();
    }

    if (tid == 0) {
        // set last element to zero for down-sweep phase
        swap(s_array[boff(2 * threads - 1)], block_sum);
    }
    __syncthreads();

    // down-sweep phase from root to leaves of binary tree
    for (uint d = 1, n = threads; n > 0; d <<= 1, n >>= 1) {
        if (tid < d) {
            uint const i1 = boff(n * (2 * tid + 1) - 1);
            uint const i2 = boff(n * (2 * tid + 2) - 1);
            T const t1 = s_array[i1];
            T const t2 = s_array[i2];
            s_array[i1] = t2;
            s_array[i2] = t1 + t2;
        }
        __syncthreads();
    }

    // write partial prefix sums to global memory
    if (i1 < count)
        g_out[i1] = s_array[boff(tid)];
    if (i2 < count)
        g_out[i2] = s_array[boff(threads + tid)];

    // block sum for last thread in block, otherwise zero
    return block_sum;
}

/**
 * blockwise parallel exclusive prefix sum
 */
template <typename T>
__global__ void grid_prefix_sum(T const* g_in, T* g_out, T* g_block_sum, uint const count)
{
    uint const tid = threadIdx.x;
    uint const bid = blockIdx.x;

    T const block_sum =  grid_prefix_sum(g_in, g_out, count);

    if (tid == 0) {
        g_block_sum[bid] = block_sum;
    }
}

/**
 * single-block parallel exclusive prefix sum
 */
template <typename T>
__global__ void block_prefix_sum(T const* g_in, T* g_out, uint const count)
{
    grid_prefix_sum(g_in, g_out, count);
}

/**
 * add block prefix sum to partial prefix sums for each block
 */
template <typename T>
__global__ void add_block_sums(T const* g_in, T* g_out, T const* g_block_sum, uint const count)
{
    __shared__ T s_block_sum[1];

    uint const tid = threadIdx.x;
    uint const threads = blockDim.x;
    uint const bid = blockIdx.x;

    if (tid == 0) {
        // read block sum for subsequent shared memory broadcast
        s_block_sum[0] = g_block_sum[bid];
    }
    __syncthreads();

    uint const i1 = 2 * bid * threads + tid;
    if (i1 < count)
        g_out[i1] = g_in[i1] + s_block_sum[0];

    uint const i2 = (2 * bid + 1) * threads + tid;
    if (i2 < count)
        g_out[i2] = g_in[i2] + s_block_sum[0];
}

} // namespace scan_kernel

/**
 * CUDA C++ wrapper
 */
template <typename T>
scan_wrapper<T> const scan_wrapper<T>::kernel = {
    scan_kernel::grid_prefix_sum
  , scan_kernel::add_block_sums
  , scan_kernel::block_prefix_sum
};

}} // namespace algorithm::gpu

} // namespace halmd

#endif /* ! HALMD_ALGORITHM_GPU_SCAN_KERNEL_CUH */
