/*
 * Copyright © 2021       Jaslo Ziska
 * Copyright © 2019       Felix Höfling
 * Copyright © 2008-2009  Peter Colberg
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

#ifndef HALMD_ALGORITHM_GPU_REDUCTION_CUH
#define HALMD_ALGORITHM_GPU_REDUCTION_CUH

#include <halmd/algorithm/gpu/bits/shfl.cuh>
#include <halmd/algorithm/gpu/transform.cuh>
#include <halmd/utility/gpu/thread.cuh>

namespace halmd {
namespace algorithm {
namespace gpu {

using halmd::algorithm::gpu::bits::shfl_down;

/**
 * parallel unary warp reduction
 */
template <typename transform_, typename T>
__device__ void reduce_warp(T& val)
{
    int delta = WARP_SIZE / 2;
    val = transform<transform_>(val, shfl_down(FULL_MASK, val, delta));

    if (WTID < WARP_SIZE / 2) { // exit upper half-warp
        while (delta > 1) {
            delta >>= 1;  // <=> delta /= 2;
            val = transform<transform_>(val, shfl_down(0xFFFF, val, delta));
        }
    }
}

/**
 * parallel unary warp reduction for less than 32 threads
 */
template <typename transform_, typename T>
__device__ void reduce_warp(T& val, unsigned int size)
{
    assert(size > 1 && size <= warpSize);

    int delta = 1 << (31 - __clz(size - 1)); // round down to previous power of 2 (powers of two are rounded down too)

    T tmp = shfl_down(FULL_MASK, val, delta);
    if (WTID + delta < size) { // only apply transformation if the warps are included in size
        val = transform<transform_>(val, tmp);
    }

    if (WTID < (warpSize >> 1)) { // exit upper half-warp
        // now that delta is a power of two reduce normally
        while (delta > 1) {
            delta >>= 1; // <=> delta /= 2;
            val = transform<transform_>(val, shfl_down(0xFFFF, val, delta));
        }
    }
}

/**
 * parallel unary reduction
 */
template <typename transform_, typename T>
__device__ void reduce(T& val)
{
    /*
     * The size of the shared memory array is 1024 / WARP_SIZE because, at the time of writing this function,
     * the maximum number of threads per block is 1024 and the warp size is 32 so the length of the shared memory
     * would never be exceeded.
     */
    assert(WDIM <= 1024 / WARP_SIZE);
    static __shared__ unsigned char shared_bytes[sizeof(T) * 1024 / WARP_SIZE];
    T* const shared = reinterpret_cast<T*>(shared_bytes);

    // reduce each warp individually
    reduce_warp<transform_>(val);

    if (WDIM > 1) { // only reduce the results if more than one warp was reduced
        // write the result from the first thread of each warp (first lane) to shared memory
        if (WTID == 0) {
            shared[WID] = val;
        }

        __syncthreads();

        // reduce the results from all warps in the first warp
        if (WID == 0) {
            val = shared[TID];
            reduce_warp<transform_>(val, WDIM);
        }
    }
}

/**
 * parallel binary warp reduction
 */
template <typename transform_, typename T0, typename T1>
__device__ void reduce_warp(T0& val0, T1& val1)
{
    int delta = WARP_SIZE / 2;
    transform<transform_>(val0, val1, shfl_down(FULL_MASK, val0, delta), shfl_down(FULL_MASK, val1, delta));

    if (WTID < WARP_SIZE / 2) { // exit upper half-warp
        while (delta > 1) {
            delta >>= 1;  // <=> delta /= 2;
            transform<transform_>(val0, val1, shfl_down(0xFFFF, val0, delta), shfl_down(0xFFFF, val1, delta));
        }
    }
}

/**
 * parallel binary warp reduction for less than 32 threads
 */
template <typename transform_, typename T0, typename T1>
__device__ void reduce_warp(T0& val0, T1& val1, unsigned int size)
{
    assert(size > 1 && size <= warpSize);

    int delta = 1 << (31 - __clz(size - 1)); // round down to previous power of 2 (powers of two are rounded down too)

    T0 tmp0 = shfl_down(FULL_MASK, val0, delta);
    T1 tmp1 = shfl_down(FULL_MASK, val1, delta);
    if (WTID + delta < size) { // only apply transformation if the warps are included in size
        transform<transform_>(val0, val1, tmp0, tmp1);
    }

    if (WTID < (warpSize >> 1)) { // exit upper half-warp
        // now that delta is a power of two reduce normally
        while (delta > 1) {
            delta = delta >> 1;
            transform<transform_>(val0, val1, shfl_down(0xFFFF, val0, delta), shfl_down(0xFFFF, val1, delta));
        }
    }
}

/**
 * parallel binary reduction
 */
template <typename transform_, typename T0, typename T1>
__device__ void reduce(T0& val0, T1& val1)
{
    /*
     * The size of the shared memory array is 1024 / WARP_SIZE because, at the time of writing this function,
     * the maximum number of threads per block is 1024 and the warp size is 32 so the length of the shared memory
     * would never be exceeded.
     */
    assert(WDIM <= 1024 / WARP_SIZE);
    static __shared__ unsigned char shared_bytes[(sizeof(T0) + sizeof(T1)) * 1024 / WARP_SIZE];
    T0* const shared0 = reinterpret_cast<T0*>(shared_bytes);
    T1* const shared1 = reinterpret_cast<T1*>(&shared0[WARP_SIZE]);

    // reduce each warp individually
    reduce_warp<transform_>(val0, val1);

    if (WDIM > 1) { // only reduce the results if more than one warp was reduced
        // write the results from the first thread of each warp (first lane) to shared memory
        if (WTID == 0) {
            shared0[WID] = val0;
            shared1[WID] = val1;
        }

        __syncthreads();

        // reduce the results from all warps in the first warp
        if (WID == 0) {
            val0 = shared0[TID];
            val1 = shared1[TID];
            reduce_warp<transform_>(val0, val1);
        }
    }
}

/**
 * parallel ternary warp reduction
 */
template <typename transform_, typename T0, typename T1, typename T2>
__device__ void reduce_warp(T0& val0, T1& val1, T2& val2)
{
    int delta = WARP_SIZE / 2;
    transform<transform_>(val0, val1, val2
      , shfl_down(FULL_MASK, val0, delta)
      , shfl_down(FULL_MASK, val1, delta)
      , shfl_down(FULL_MASK, val2, delta)
    );

    if (WTID < WARP_SIZE / 2) { // exit upper half-warp
        while (delta > 1) {
            delta >>= 1;  // <=> delta /= 2;
            transform<transform_>(val0, val1, val2
              , shfl_down(0xFFFF, val0, delta)
              , shfl_down(0xFFFF, val1, delta)
              , shfl_down(0xFFFF, val2, delta)
            );
        }
    }
}

/**
 * parallel ternary warp reduction for less than 32 threads
 */
template <typename transform_, typename T0, typename T1, typename T2>
__device__ void reduce_warp(T0& val0, T1& val1, T2& val2, unsigned int size)
{
    assert(size > 1 && size <= warpSize);

    int delta = 1 << (31 - __clz(size - 1)); // round down to previous power of 2 (powers of two are rounded down too)

    T0 tmp0 = shfl_down(FULL_MASK, val0, delta);
    T1 tmp1 = shfl_down(FULL_MASK, val1, delta);
    T2 tmp2 = shfl_down(FULL_MASK, val2, delta);
    if (WTID + delta < size) { // only apply transformation if the warps are included in size
        transform<transform_>(val0, val1, val2, tmp0, tmp1, tmp2);
    }

    if (WTID < (warpSize >> 1)) { // exit upper half-warp
        // now that delta is a power of two reduce normally
        while (delta > 1) {
            delta = delta >> 1;
            transform<transform_>(val0, val1, val2
              , shfl_down(0xFFFF, val0, delta)
              , shfl_down(0xFFFF, val1, delta)
              , shfl_down(0xFFFF, val2, delta)
            );
        }
    }
}


/**
 * parallel ternary reduction
 */
template <typename transform_, typename T0, typename T1, typename T2>
__device__ void reduce(T0& val0, T1& val1, T2& val2)
{
    /*
     * The size of the shared memory array is 1024 / WARP_SIZE because, at the time of writing this function,
     * the maximum number of threads per block is 1024 and the warp size is 32 so the length of the shared memory
     * would never be exceeded.
     */
    assert(WDIM <= 1024 / WARP_SIZE);
    static __shared__ unsigned char shared_bytes[(sizeof(T0) + sizeof(T1) + sizeof(T2)) * 1024 / WARP_SIZE];
    T0* const shared0 = reinterpret_cast<T0*>(shared_bytes);
    T1* const shared1 = reinterpret_cast<T1*>(&shared0[WARP_SIZE]);
    T2* const shared2 = reinterpret_cast<T2*>(&shared1[WARP_SIZE]);

    // reduce each warp individually
    reduce_warp<transform_>(val0, val1, val2);

    if (WDIM > 1) { // only reduce the results if more than one warp was reduced
        // write the results from the first thread of each warp (first lane) to shared memory
        if (WTID == 0) {
            shared0[WID] = val0;
            shared1[WID] = val1;
            shared2[WID] = val2;
        }

        __syncthreads();

        // reduce the results from all warps in the first warp
        if (WID == 0) {
            val0 = shared0[WTID];
            val1 = shared1[WTID];
            val2 = shared2[WTID];
            reduce_warp<transform_>(val0, val1, val2, WDIM);
        }
    }
}

/**
 * parallel quartenary warp reduction
 */
template <typename transform_, typename T0, typename T1, typename T2, typename T3>
__device__ void reduce_warp(T0& val0, T1& val1, T2& val2, T3& val3)
{
    int delta = WARP_SIZE / 2;
    transform<transform_>(val0, val1, val2, val3
      , shfl_down(FULL_MASK, val0, delta)
      , shfl_down(FULL_MASK, val1, delta)
      , shfl_down(FULL_MASK, val2, delta)
      , shfl_down(FULL_MASK, val3, delta)
    );

    if (WTID < WARP_SIZE / 2) { // exit upper half-warp
        while (delta > 1) {
            delta >>= 1;  // <=> delta /= 2;
            transform<transform_>(val0, val1, val2, val3
              , shfl_down(0xFFFF, val0, delta)
              , shfl_down(0xFFFF, val1, delta)
              , shfl_down(0xFFFF, val2, delta)
              , shfl_down(0xFFFF, val3, delta)
            );
        }
    }
}

/**
 * parallel quartenary warp reduction for less than 32 threads
 */
template <typename transform_, typename T0, typename T1, typename T2, typename T3>
__device__ void reduce_warp(T0& val0, T1& val1, T2& val2, T3& val3, unsigned int size)
{
    assert(size > 1 && size <= warpSize);

    int delta = 1 << (31 - __clz(size - 1)); // round down to previous power of 2 (powers of two are rounded down too)

    T0 tmp0 = shfl_down(FULL_MASK, val0, delta);
    T1 tmp1 = shfl_down(FULL_MASK, val1, delta);
    T2 tmp2 = shfl_down(FULL_MASK, val2, delta);
    T3 tmp3 = shfl_down(FULL_MASK, val3, delta);
    if (WTID + delta < size) { // only apply transformation if the warps are included in size
        transform<transform_>(val0, val1, val2, val3, tmp0, tmp1, tmp2, tmp3);
    }

    if (WTID < (warpSize >> 1)) { // exit upper half-warp
        // now that delta is a power of two reduce normally
        while (delta > 1) {
            delta = delta >> 1;
            transform<transform_>(val0, val1, val2, val3
              , shfl_down(0xFFFF, val0, delta)
              , shfl_down(0xFFFF, val1, delta)
              , shfl_down(0xFFFF, val2, delta)
              , shfl_down(0xFFFF, val3, delta)
            );
        }
    }
}


/**
 * parallel quartenary reduction
 */
template <typename transform_, typename T0, typename T1, typename T2, typename T3>
__device__ void reduce(T0& val0, T1& val1, T2& val2, T3& val3)
{
    /*
     * The size of the shared memory array is 1024 / WARP_SIZE because, at the time of writing this function,
     * the maximum number of threads per block is 1024 and the warp size is 32 so the length of the shared memory
     * would never be exceeded.
     */
    assert(WDIM <= 1024 / WARP_SIZE);
    static __shared__ unsigned char shared_bytes[(sizeof(T0) + sizeof(T1) + sizeof(T2) + sizeof(T3)) * 1024 / WARP_SIZE];
    T0* const shared0 = reinterpret_cast<T0*>(shared_bytes);
    T1* const shared1 = reinterpret_cast<T1*>(&shared0[WARP_SIZE]);
    T2* const shared2 = reinterpret_cast<T2*>(&shared1[WARP_SIZE]);
    T3* const shared3 = reinterpret_cast<T3*>(&shared2[WARP_SIZE]);

    // reduce each warp individually
    reduce_warp<transform_>(val0, val1, val2, val3);

    if (WDIM > 1) { // only reduce the results if more than one warp was reduced
        // write the results from the first thread of each warp (first lane) to shared memory
        if (WTID == 0) {
            shared0[WID] = val0;
            shared1[WID] = val1;
            shared2[WID] = val2;
            shared3[WID] = val3;
        }

        __syncthreads();

        // reduce the results from all warps in the first warp
        if (WID == 0) {
            val0 = shared0[TID];
            val1 = shared1[TID];
            val2 = shared2[TID];
            val3 = shared3[TID];
            reduce_warp<transform_>(val0, val1, val2, val3, WDIM);
        }
    }
}

} // namespace algorithm
} // namespace gpu
} // namespace halmd

#endif /* ! HALMD_ALGORITHM_GPU_REDUCTION_CUH */
