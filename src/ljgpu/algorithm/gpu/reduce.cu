/* Parallel reduction kernel
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

#include <ljgpu/algorithm/gpu/base.cuh>
#include <ljgpu/algorithm/gpu/reduce.hpp>
#include <ljgpu/math/gpu/vector2d.cuh>
#include <ljgpu/math/gpu/vector3d.cuh>

namespace ljgpu { namespace cu { namespace reduce
{

enum { THREADS = gpu::reduce::THREADS };

/**
 * unary transformations
 */
template <typename T, typename U>
__device__ T identity_(U v)
{
    return v;
}

template <>
__device__ float3 identity_(float4 v)
{
    return make_float3(v.x, v.y, v.z);
}

template <>
__device__ float4 identity_(float3 v)
{
    return make_float4(v.x, v.y, v.z, 0);
}

template <typename T, typename U>
__device__ T square_(U v)
{
    return v * v;
}

/**
 * binary transformations
 */
template <typename T>
__device__ T sum_(T v1, T v2)
{
    return v1 + v2;
}

/**
 * parallel reduction
 */
template <typename input_type, typename output_type,
	  output_type (*input_function)(input_type),
	  output_type (*reduce_function)(output_type, output_type),
	  output_type (*output_function)(output_type),
	  typename T, typename U>
__device__ void reduce(T const* g_in, U* g_block_sum, uint n)
{
    __shared__ output_type s_vv[THREADS];

    // load values from global device memory
    output_type vv = 0;
    for (uint i = GTID; i < n; i += GTDIM) {
	input_type v = identity_<input_type>(g_in[i]);
	vv = reduce_function(vv, input_function(v));
    }
    // reduced value for this thread
    s_vv[TID] = vv;
    __syncthreads();

    // compute reduced value for all threads in block
    if (TID < 256) {
	vv = reduce_function(vv, s_vv[TID + 256]);
	s_vv[TID] = vv;
    }
    __syncthreads();
    if (TID < 128) {
	vv = reduce_function(vv, s_vv[TID + 128]);
	s_vv[TID] = vv;
    }
    __syncthreads();
    if (TID < 64) {
	vv = reduce_function(vv, s_vv[TID + 64]);
	s_vv[TID] = vv;
    }
    __syncthreads();
    if (TID < 32) {
	vv = reduce_function(vv, s_vv[TID + 32]);
	s_vv[TID] = vv;
    }
    // no further syncs needed within execution warp of 32 threads
    if (TID < 16) {
	vv = reduce_function(vv, s_vv[TID + 16]);
	s_vv[TID] = vv;
    }
    if (TID < 8) {
	vv = reduce_function(vv, s_vv[TID + 8]);
	s_vv[TID] = vv;
    }
    if (TID < 4) {
	vv = reduce_function(vv, s_vv[TID + 4]);
	s_vv[TID] = vv;
    }
    if (TID < 2) {
	vv = reduce_function(vv, s_vv[TID + 2]);
	s_vv[TID] = vv;
    }
    if (TID < 1) {
	vv = reduce_function(vv, s_vv[TID + 1]);
	// store block reduced value in global memory
	g_block_sum[blockIdx.x] = identity_<U>(output_function(vv));
    }
}

/**
 * blockwise sum
 */
template <typename input_type, typename output_type, typename T, typename U>
__global__ void sum(T const* g_in, U* g_block_sum, uint n)
{
    reduce<input_type, output_type, identity_, sum_, identity_>(g_in, g_block_sum, n);
}

/**
 * blockwise absolute maximum
 */
template <typename input_type, typename output_type, typename T, typename U>
__global__ void max(T const* g_in, U* g_block_max, uint n)
{
    reduce<input_type, output_type, square_, fmaxf, sqrtf>(g_in, g_block_max, n);
}

}}} // namespace ljgpu::cu::reduce

namespace ljgpu { namespace gpu
{

/**
 * device function wrappers
 */
cuda::function<void(float const*, dfloat*, uint)>
	       reduce::sum(cu::reduce::sum<float, dfloat>);
cuda::function<void(float4 const*, float*, uint),
	       void(float2 const*, float*, uint)>
	       reduce::max(cu::reduce::max<float3, float>,
			   cu::reduce::max<float2, float>);

}} // namespace ljgpu::gpu
