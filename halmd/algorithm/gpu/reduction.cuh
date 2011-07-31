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

#ifndef HALMD_ALGORITHM_GPU_REDUCTION_CUH
#define HALMD_ALGORITHM_GPU_REDUCTION_CUH

#include <boost/mpl/int.hpp>

#include <halmd/algorithm/gpu/transform.cuh>

namespace halmd {
namespace algorithm {
namespace gpu {

using boost::mpl::int_;

/**
 * parallel unary reduction
 */
template <int threads, typename transform_, typename T, typename V>
__device__ typename enable_if<is_same<int_<threads>, int_<1> >, void>::type
reduce(T& sum, V s_sum[])
{
    int const tid = threadIdx.x;
    if (tid < threads) {
        sum = transform<transform_>(sum, static_cast<T>(s_sum[tid + threads]));
    }
}

template <int threads, typename transform_, typename T, typename V>
__device__ typename disable_if<is_same<int_<threads>, int_<1> >, void>::type
reduce(T& sum, V s_sum[])
{
    int const tid = threadIdx.x;
    if (tid < threads) {
        sum = transform<transform_>(sum, static_cast<T>(s_sum[tid + threads]));
        s_sum[tid] = sum;
    }
    // no further syncs needed within execution warp of 32 threads
    if (threads >= warpSize) {
        __syncthreads();
    }

    reduce<threads / 2, transform_>(sum, s_sum);
}

/**
 * parallel binary reduction
 */
template <int threads, typename transform_, typename T0, typename T1, typename V0, typename V1>
__device__ typename enable_if<is_same<int_<threads>, int_<1> >, void>::type
reduce(T0& sum0, T1& sum1, V0 s_sum0[], V1 s_sum1[])
{
    int const tid = threadIdx.x;
    if (tid < threads) {
        transform<transform_>(sum0, sum1, static_cast<T0>(s_sum0[tid + threads]), static_cast<T1>(s_sum1[tid + threads]));
    }
}

template <int threads, typename transform_, typename T0, typename T1, typename V0, typename V1>
__device__ typename disable_if<is_same<int_<threads>, int_<1> >, void>::type
reduce(T0& sum0, T1& sum1, V0 s_sum0[], V1 s_sum1[])
{
    int const tid = threadIdx.x;
    if (tid < threads) {
        transform<transform_>(sum0, sum1, static_cast<T0>(s_sum0[tid + threads]), static_cast<T1>(s_sum1[tid + threads]));
        s_sum0[tid] = sum0;
        s_sum1[tid] = sum1;
    }
    // no further syncs needed within execution warp of 32 threads
    if (threads >= warpSize) {
        __syncthreads();
    }

    reduce<threads / 2, transform_>(sum0, sum1, s_sum0, s_sum1);
}

/**
 * parallel ternary reduction
 */
template <int threads, typename transform_, typename T0, typename T1, typename T2, typename V0, typename V1, typename V2>
__device__ typename enable_if<is_same<int_<threads>, int_<1> >, void>::type
reduce(T0& sum0, T1& sum1, T2& sum2, V0 s_sum0[], V1 s_sum1[], V2 s_sum2[])
{
    int const tid = threadIdx.x;
    if (tid < threads) {
        transform<transform_>(sum0, sum1, sum2, static_cast<T0>(s_sum0[tid + threads]), static_cast<T1>(s_sum1[tid + threads]), static_cast<T2>(s_sum2[tid + threads]));
    }
}

template <int threads, typename transform_, typename T0, typename T1, typename T2, typename V0, typename V1, typename V2>
__device__ typename disable_if<is_same<int_<threads>, int_<1> >, void>::type
reduce(T0& sum0, T1& sum1, T2& sum2, V0 s_sum0[], V1 s_sum1[], V2 s_sum2[])
{
    int const tid = threadIdx.x;
    if (tid < threads) {
        transform<transform_>(sum0, sum1, sum2, static_cast<T0>(s_sum0[tid + threads]), static_cast<T1>(s_sum1[tid + threads]), static_cast<T2>(s_sum2[tid + threads]));
        s_sum0[tid] = sum0;
        s_sum1[tid] = sum1;
        s_sum2[tid] = sum2;
    }
    // no further syncs needed within execution warp of 32 threads
    if (threads >= warpSize) {
        __syncthreads();
    }

    reduce<threads / 2, transform_>(sum0, sum1, sum2, s_sum0, s_sum1, s_sum2);
}

/**
 * parallel quartenary reduction
 */
template <int threads, typename transform_, typename T0, typename T1, typename T2, typename T3, typename V0, typename V1, typename V2, typename V3>
__device__ typename enable_if<is_same<int_<threads>, int_<1> >, void>::type
reduce(T0& sum0, T1& sum1, T2& sum2, T3& sum3, V0 s_sum0[], V1 s_sum1[], V2 s_sum2[], V3 s_sum3[])
{
    int const tid = threadIdx.x;
    if (tid < threads) {
        transform<transform_>(sum0, sum1, sum2, sum3, static_cast<T0>(s_sum0[tid + threads]), static_cast<T1>(s_sum1[tid + threads]), static_cast<T2>(s_sum2[tid + threads]), static_cast<T3>(s_sum3[tid + threads]));
    }
}

template <int threads, typename transform_, typename T0, typename T1, typename T2, typename T3, typename V0, typename V1, typename V2, typename V3>
__device__ typename disable_if<is_same<int_<threads>, int_<1> >, void>::type
reduce(T0& sum0, T1& sum1, T2& sum2, T3& sum3, V0 s_sum0[], V1 s_sum1[], V2 s_sum2[], V3 s_sum3[])
{
    int const tid = threadIdx.x;
    if (tid < threads) {
        transform<transform_>(sum0, sum1, sum2, sum3, static_cast<T0>(s_sum0[tid + threads]), static_cast<T1>(s_sum1[tid + threads]), static_cast<T2>(s_sum2[tid + threads]), static_cast<T3>(s_sum3[tid + threads]));
        s_sum0[tid] = sum0;
        s_sum1[tid] = sum1;
        s_sum2[tid] = sum2;
        s_sum3[tid] = sum3;
    }
    // no further syncs needed within execution warp of 32 threads
    if (threads >= warpSize) {
        __syncthreads();
    }

    reduce<threads / 2, transform_>(sum0, sum1, sum2, sum3, s_sum0, s_sum1, s_sum2, s_sum3);
}

} // namespace algorithm
} // namespace gpu
} // namespace halmd

#endif /* ! HALMD_ALGORITHM_GPU_REDUCTION_CUH */
