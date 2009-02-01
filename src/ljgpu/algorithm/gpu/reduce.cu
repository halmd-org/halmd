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
#include <ljgpu/algorithm/gpu/reduce.cuh>
#include <ljgpu/algorithm/gpu/reduce.hpp>
#include <ljgpu/math/gpu/vector2d.cuh>
#include <ljgpu/math/gpu/vector3d.cuh>

namespace ljgpu { namespace cu { namespace algorithm
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
	  typename coalesced_input_type, typename coalesced_output_type>
__device__ void accumulate(coalesced_input_type const* g_in, coalesced_output_type* g_block_sum, uint n)
{
    __shared__ output_type s_vv[THREADS];

    // load values from global device memory
    output_type vv = 0;
    for (uint i = GTID; i < n; i += GTDIM) {
	input_type v = g_in[i];
	vv = reduce_function(vv, input_function(v));
    }
    // reduced value for this thread
    s_vv[TID] = vv;
    __syncthreads();

    // compute reduced value for all threads in block
    reduce<THREADS / 2>(vv, s_vv);

    if (TID < 1) {
	// store block reduced value in global memory
	g_block_sum[blockIdx.x] = output_function(vv);
    }
}

/**
 * blockwise sum
 */
template <typename input_type, typename output_type,
	  typename coalesced_input_type, typename coalesced_output_type>
__global__ void sum(coalesced_input_type const* g_in, coalesced_output_type* g_block_sum, uint n)
{
    accumulate<input_type, output_type, identity_, sum_, identity_>(g_in, g_block_sum, n);
}

/**
 * blockwise sum of squares
 */
template <typename input_type, typename output_type,
	  typename coalesced_input_type, typename coalesced_output_type>
__global__ void sum_of_squares(coalesced_input_type const* g_in, coalesced_output_type* g_block_sum, uint n)
{
    accumulate<input_type, output_type, square_, sum_, identity_>(g_in, g_block_sum, n);
}

/**
 * blockwise absolute maximum
 */
template <typename input_type, typename output_type,
	  typename coalesced_input_type, typename coalesced_output_type>
__global__ void max(coalesced_input_type const* g_in, coalesced_output_type* g_block_max, uint n)
{
    accumulate<input_type, output_type, square_, fmaxf, sqrtf>(g_in, g_block_max, n);
}

}}} // namespace ljgpu::cu::algorithm

namespace ljgpu { namespace gpu
{

/**
 * device function wrappers
 */
cuda::function<void(float const*, dfloat*, uint),
	       void(float4 const*, float4*, uint),
	       void(float2 const*, float2*, uint)>
    reduce::sum(cu::algorithm::sum<float, dfloat>,
		cu::algorithm::sum<cu::vector<float, 3>, cu::vector<float, 3> >,
		cu::algorithm::sum<cu::vector<float, 2>, cu::vector<float, 2> >);
cuda::function<void(float4 const*, dfloat*, uint),
	       void(float2 const*, dfloat*, uint)>
    reduce::sum_of_squares(cu::algorithm::sum_of_squares<cu::vector<float, 3>, dfloat>,
			   cu::algorithm::sum_of_squares<cu::vector<float, 2>, dfloat>);
cuda::function<void(float4 const*, float*, uint),
	       void(float2 const*, float*, uint)>
    reduce::max(cu::algorithm::max<cu::vector<float, 3>, float>,
		cu::algorithm::max<cu::vector<float, 2>, float>);

}} // namespace ljgpu::gpu
