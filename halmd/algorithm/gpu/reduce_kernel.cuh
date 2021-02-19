/*
 * Copyright © 2008-2012 Peter Colberg
 * Copyright © 2021      Jaslo Ziska
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

#ifndef HALMD_ALGORITHM_GPU_REDUCE_KERNEL_CUH
#define HALMD_ALGORITHM_GPU_REDUCE_KERNEL_CUH

#include <boost/utility/enable_if.hpp>

#include <halmd/algorithm/gpu/bits/shfl.cuh>
#include <halmd/algorithm/gpu/reduce_kernel.hpp>
#include <halmd/utility/gpu/thread.cuh>

namespace halmd {
namespace detail {

using halmd::algorithm::gpu::bits::shfl_down;

/**
 * Reduce accumulators of warp
 *
 * @param acc accumulator
 */
template <typename accumulator_type>
__device__ void reduce_warp(accumulator_type& acc)
{
    int delta = WARP_SIZE / 2;
    acc(shfl_down(FULL_MASK, acc, delta));

    if (WTID < WARP_SIZE / 2) {
        while (delta > 1) {
            delta >>= 1; // <=> delta /= 2;
            acc(shfl_down(0xFFFF, acc, delta));
        }
    }
}

/**
 * Reduce accumulators of warp for less than 32 threads
 *
 * @param acc accumulator
 * @param size number of threads to be reduced
 */
template <typename accumulator_type>
__device__ void reduce_warp(accumulator_type& acc, unsigned int size)
{
    assert(size > 1 && size <= warpSize);

    int delta = 1 << (31 - __clz(size - 1)); // round down to previous power of 2 (powers of two are rounded down too)

    accumulator_type tmp = shfl_down(FULL_MASK, acc, delta);
    if (WTID + delta < size) { // only apply transformation if the warps are included in size
        acc(tmp);
    }

    if (WTID < (warpSize >> 1)) { // exit upper half-warp
        // now that delta is a power of two reduce normally
        while (delta > 1) {
            delta >>= 1; // <=> delta /= 2;
            acc(shfl_down(0xFFFF, acc, delta));
        }
    }
}

/**
 * Reduce accumulators of block threads.
 *
 * @param acc accumulator of block thread 0
 */
template <typename accumulator_type>
inline __device__ void reduce(accumulator_type& acc)
{
    // We need to avoid default initialization of the shared memory
    // array, since this increases execution time of the kernel.
    // Use a char array, and cast to type of reduction functor.
    // The size of the shared memory array is 1024 / WARP_SIZE because, at the time of writing this function,
    // the maximum number of threads per block is 1024 and the warp size is 32 so the length of the shared memory
    // would never be exceeded.
    assert(WDIM <= 1024 / WARP_SIZE);
    __shared__ char s_storage[sizeof(accumulator_type) * 1024 / WARP_SIZE];
    accumulator_type* const s_acc = reinterpret_cast<accumulator_type*>(s_storage);

    // reduce each warp individually
    reduce_warp(acc);

    // write the result from the first thread of each warp (first lane) to shared memory
    if (WTID == 0) {
        s_acc[WID] = acc;
    }

    __syncthreads();

    // reduce the results from all warps in the first warp
    if (WID == 0) {
        acc = s_acc[WTID];
        reduce_warp(acc, WDIM);
    }
}

/**
 * Compute block sums of input array using unary accumulator.
 *
 * @param first input iterator to first element
 * @param size number of elements
 * @param g_block_acc output block accumulators
 * @param input accumulator
 */
template <typename accumulator_type>
static __global__ void reduction(
    typename accumulator_type::iterator const first
  , typename std::iterator_traits<typename accumulator_type::iterator>::difference_type size
  , accumulator_type* g_block_acc
  , accumulator_type acc
)
{
    // load values from global device memory
    for (unsigned int i = GTID; i < size; i += GTDIM) {
        acc(first[i]);
    }
    // compute reduced value for all threads in block
    reduce(acc);

    if (TID < 1) {
        // store block reduced value in global memory
        g_block_acc[blockIdx.x] = acc;
    }
}

} // namespace detail

template <typename accumulator_type>
reduction_kernel<accumulator_type> reduction_kernel<accumulator_type>::kernel = {
    detail::reduction
};

} // namespace halmd

#endif /* ! HALMD_ALGORITHM_GPU_REDUCE_KERNEL_CUH */
