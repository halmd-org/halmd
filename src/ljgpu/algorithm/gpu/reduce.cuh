/* Parallel reduction
 *
 * Copyright © 2008-2009  Peter Colberg
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

#ifndef LJGPU_ALGORITHM_GPU_REDUCE_CUH
#define LJGPU_ALGORITHM_GPU_REDUCE_CUH

#include <boost/mpl/int.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>
#include <ljgpu/math/gpu/accum.cuh>
using namespace boost;

namespace ljgpu { namespace cu
{

enum { WARP_SIZE = 32 };

/**
 * transformation types
 */
struct identity_;
struct square_;
struct sqrt_;
struct sum_;
struct complex_sum_;
struct quaternion_sum_;
struct max_;

/**
 * unary transformations
 */
template <typename transform_, typename input_type, typename output_type>
__device__ typename enable_if<is_same<transform_, identity_>, output_type>::type
transform(input_type v)
{
    return v;
}

template <typename transform_, typename input_type, typename output_type>
__device__ typename enable_if<is_same<transform_, square_>, output_type>::type
transform(input_type v)
{
    return v * v;
}

template <typename transform_, typename input_type, typename output_type>
__device__ typename enable_if<is_same<transform_, sqrt_>, output_type>::type
transform(input_type v)
{
    return sqrtf(v);
}

/**
 * binary transformations
 */
template <typename transform_, typename T>
__device__ typename enable_if<is_same<transform_, sum_>, T>::type
transform(T v1, T v2)
{
    return v1 + v2;
}

template <typename transform_, typename T>
__device__ typename enable_if<is_same<transform_, complex_sum_>, void>::type
transform(T& r1, T& i1, T r2, T i2)
{
    r1 += r2;
    i1 += i2;
}

template <typename transform_, typename T>
__device__ typename enable_if<is_same<transform_, quaternion_sum_>, void>::type
transform(T& r1, T& i1, T& j1, T& k1, T r2, T i2, T j2, T k2)
{
    r1 += r2;
    i1 += i2;
    j1 += j2;
    k1 += k2;
}

template <typename transform_, typename T>
__device__ typename enable_if<is_same<transform_, max_>, T>::type
transform(T v1, T v2)
{
    return fmaxf(v1, v2);
}

/**
 * parallel unary reduction
 */
template <int threads, typename transform_, typename T>
__device__ typename enable_if<is_same<mpl::int_<threads>, mpl::int_<1> >, void>::type
reduce(T& sum, T s_sum[])
{
    int const tid = threadIdx.x;
    if (tid < threads) {
	sum = transform<transform_>(sum, s_sum[tid + threads]);
    }
}

template <int threads, typename transform_, typename T>
__device__ typename disable_if<is_same<mpl::int_<threads>, mpl::int_<1> >, void>::type
reduce(T& sum, T s_sum[])
{
    int const tid = threadIdx.x;
    if (tid < threads) {
	sum = transform<transform_>(sum, s_sum[tid + threads]);
	s_sum[tid] = sum;
    }
    // no further syncs needed within execution warp of 32 threads
    if (threads >= WARP_SIZE) {
	__syncthreads();
    }

    reduce<threads / 2, transform_>(sum, s_sum);
}

/**
 * parallel binary reduction
 */
template <int threads, typename transform_, typename T0, typename T1>
__device__ typename enable_if<is_same<mpl::int_<threads>, mpl::int_<1> >, void>::type
reduce(T0& sum0, T1& sum1, T0 s_sum0[], T1 s_sum1[])
{
    int const tid = threadIdx.x;
    if (tid < threads) {
	transform<transform_>(sum0, sum1, s_sum0[tid + threads], s_sum1[tid + threads]);
    }
}

template <int threads, typename transform_, typename T0, typename T1>
__device__ typename disable_if<is_same<mpl::int_<threads>, mpl::int_<1> >, void>::type
reduce(T0& sum0, T1& sum1, T0 s_sum0[], T1 s_sum1[])
{
    int const tid = threadIdx.x;
    if (tid < threads) {
	transform<transform_>(sum0, sum1, s_sum0[tid + threads], s_sum1[tid + threads]);
	s_sum0[tid] = sum0;
	s_sum1[tid] = sum1;
    }
    // no further syncs needed within execution warp of 32 threads
    if (threads >= WARP_SIZE) {
	__syncthreads();
    }

    reduce<threads / 2, transform_>(sum0, sum1, s_sum0, s_sum1);
}

/**
 * parallel ternary reduction
 */
template <int threads, typename transform_, typename T0, typename T1, typename T2>
__device__ typename enable_if<is_same<mpl::int_<threads>, mpl::int_<1> >, void>::type
reduce(T0& sum0, T1& sum1, T2& sum2, T0 s_sum0[], T1 s_sum1[], T2 s_sum2[])
{
    int const tid = threadIdx.x;
    if (tid < threads) {
	transform<transform_>(sum0, sum1, sum2, s_sum0[tid + threads], s_sum1[tid + threads], s_sum2[tid + threads]);
    }
}

template <int threads, typename transform_, typename T0, typename T1, typename T2>
__device__ typename disable_if<is_same<mpl::int_<threads>, mpl::int_<1> >, void>::type
reduce(T0& sum0, T1& sum1, T2& sum2, T0 s_sum0[], T1 s_sum1[], T2 s_sum2[])
{
    int const tid = threadIdx.x;
    if (tid < threads) {
	transform<transform_>(sum0, sum1, sum2, s_sum0[tid + threads], s_sum1[tid + threads], s_sum2[tid + threads]);
	s_sum0[tid] = sum0;
	s_sum1[tid] = sum1;
	s_sum2[tid] = sum2;
    }
    // no further syncs needed within execution warp of 32 threads
    if (threads >= WARP_SIZE) {
	__syncthreads();
    }

    reduce<threads / 2, transform_>(sum0, sum1, sum2, s_sum0, s_sum1, s_sum2);
}

/**
 * parallel quartenary reduction
 */
template <int threads, typename transform_, typename T0, typename T1, typename T2, typename T3>
__device__ typename enable_if<is_same<mpl::int_<threads>, mpl::int_<1> >, void>::type
reduce(T0& sum0, T1& sum1, T2& sum2, T3& sum3, T0 s_sum0[], T1 s_sum1[], T2 s_sum2[], T3 s_sum3[])
{
    int const tid = threadIdx.x;
    if (tid < threads) {
	transform<transform_>(sum0, sum1, sum2, sum3, s_sum0[tid + threads], s_sum1[tid + threads], s_sum2[tid + threads], s_sum3[tid + threads]);
    }
}

template <int threads, typename transform_, typename T0, typename T1, typename T2, typename T3>
__device__ typename disable_if<is_same<mpl::int_<threads>, mpl::int_<1> >, void>::type
reduce(T0& sum0, T1& sum1, T2& sum2, T3& sum3, T0 s_sum0[], T1 s_sum1[], T2 s_sum2[], T3 s_sum3[])
{
    int const tid = threadIdx.x;
    if (tid < threads) {
	transform<transform_>(sum0, sum1, sum2, sum3, s_sum0[tid + threads], s_sum1[tid + threads], s_sum2[tid + threads], s_sum3[tid + threads]);
	s_sum0[tid] = sum0;
	s_sum1[tid] = sum1;
	s_sum2[tid] = sum2;
	s_sum3[tid] = sum3;
    }
    // no further syncs needed within execution warp of 32 threads
    if (threads >= WARP_SIZE) {
	__syncthreads();
    }

    reduce<threads / 2, transform_>(sum0, sum1, sum2, sum3, s_sum0, s_sum1, s_sum2, s_sum3);
}

}} // namespace ljgpu::cu

#endif /* LJGPU_ALGORITHM_GPU_REDUCE_CUH */
