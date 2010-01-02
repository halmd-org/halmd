/* Parallel exclusive prefix sum
 *
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

#include <halmd/algorithm/gpu/base.cuh>
#include <halmd/algorithm/gpu/prefix_sum.hpp>

using namespace halmd::gpu::prefix_sum;

namespace halmd { namespace cu { namespace prefix_sum
{

/**
 * blockwise parallel exclusive prefix sum
 */
template <typename T>
__device__ T grid_prefix_sum(T const* g_in, T* g_out, const uint count)
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

    extern __shared__ T s_array[];
    T block_sum = 0;

    const uint tid = threadIdx.x;
    const uint threads = blockDim.x;
    const uint bid = blockIdx.x;

    // read elements from global memory, or pad with zero
    const uint i1 = 2 * bid * threads + tid;
    const uint i2 = (2 * bid + 1) * threads + tid;
    s_array[boff(tid)] = (i1 < count) ? g_in[i1] : 0;
    s_array[boff(threads + tid)] = (i2 < count) ? g_in[i2] : 0;
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
            const uint i1 = boff(n * (2 * tid + 1) - 1);
            const uint i2 = boff(n * (2 * tid + 2) - 1);
            const T t1 = s_array[i1];
            const T t2 = s_array[i2];
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
__global__ void grid_prefix_sum(T const* g_in, T* g_out, T* g_block_sum, const uint count)
{
    const uint tid = threadIdx.x;
    const uint bid = blockIdx.x;

    const T block_sum =  grid_prefix_sum(g_in, g_out, count);

    if (tid == 0) {
        g_block_sum[bid] = block_sum;
    }
}

/**
 * single-block parallel exclusive prefix sum
 */
template <typename T>
__global__ void block_prefix_sum(T const* g_in, T* g_out, const uint count)
{
    grid_prefix_sum(g_in, g_out, count);
}

/**
 * add block prefix sum to partial prefix sums for each block
 */
template <typename T>
__global__ void add_block_sums(T const* g_in, T* g_out, T const* g_block_sum, const uint count)
{
    __shared__ T s_block_sum[1];

    const uint tid = threadIdx.x;
    const uint threads = blockDim.x;
    const uint bid = blockIdx.x;

    if (tid == 0) {
        // read block sum for subsequent shared memory broadcast
        s_block_sum[0] = g_block_sum[bid];
    }
    __syncthreads();

    const uint i1 = 2 * bid * threads + tid;
    if (i1 < count)
        g_out[i1] = g_in[i1] + s_block_sum[0];

    const uint i2 = (2 * bid + 1) * threads + tid;
    if (i2 < count)
        g_out[i2] = g_in[i2] + s_block_sum[0];
}

}}} // namespace halmd::cu::prefix_sum

namespace halmd { namespace gpu
{

/**
 * device function wrappers
 */
cuda::function<void (uint const*, uint*, uint*, const uint),
               void (uint48 const*, uint48*, uint48*, const uint)>
               prefix_sum::grid_prefix_sum(cu::prefix_sum::grid_prefix_sum,
                                           cu::prefix_sum::grid_prefix_sum);
cuda::function<void (uint const*, uint*, uint const*, const uint),
               void (uint48 const*, uint48*, uint48 const*, const uint)>
               prefix_sum::add_block_sums(cu::prefix_sum::add_block_sums,
                                          cu::prefix_sum::add_block_sums);
cuda::function<void (uint48 const*, uint48*, const uint),
               void (uint const*, uint*, const uint)>
               prefix_sum::block_prefix_sum(cu::prefix_sum::block_prefix_sum,
                                            cu::prefix_sum::block_prefix_sum);

}} // namespace halmd::gpu
