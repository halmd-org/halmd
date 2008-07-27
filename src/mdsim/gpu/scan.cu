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

#include "scan_glue.hpp"
#include "cutil.h"

namespace mdsim
{

/**
 * blockwise parallel exclusive prefix sum
 */
template <typename T>
__global__ void block_prefix_sum(T const* g_in, T* g_out, T* g_block, const uint count)
{
    //
    // Prefix Sums and Their Applications,
    // Guy E. Blelloch, November 1990, CMU-CS-90-190
    //
    // Parallel Prefix Sum (Scan) with CUDA,
    // Mark Harris, April 2007, NVIDIA
    //

    using namespace mdsim::gpu::scan;
    extern __shared__ T s_array[];

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
	// write last block prefix sum to auxiliary array
	const uint i = boff(2 * threads - 1);
	g_block[bid] = s_array[i];
	// set last element to zero for down-sweep phase
	s_array[i] = 0;
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
}

/**
 * add block prefix sum to partial prefix sums for each block
 */
template <typename T>
__global__ void add_block_sums(T const* g_block, T const* g_in, T* g_out, const uint count)
{
    __shared__ T s_block_sum[1];

    const uint tid = threadIdx.x;
    const uint threads = blockDim.x;
    const uint bid = blockIdx.x;

    if (tid == 0) {
	// read block sum for subsequent shared memory broadcast
	s_block_sum[0] = g_block[bid];
    }
    __syncthreads();

    const uint i1 = 2 * bid * threads + tid;
    const uint i2 = (2 * bid + 1) * threads + tid;
    if (i1 < count)
	g_out[i1] = g_in[i1] + s_block_sum[0];
    if (i2 < count)
	g_out[i2] = g_in[i2] + s_block_sum[0];
}

} // namespace mdsim


namespace mdsim { namespace gpu { namespace scan
{

cuda::function<void (uint const*, uint*, uint*, const uint)> block_prefix_sum(mdsim::block_prefix_sum);
cuda::function<void (uint const*, uint const*, uint*, const uint)> add_block_sums(mdsim::add_block_sums);

}}} // namespace mdsim::gpu::scan
