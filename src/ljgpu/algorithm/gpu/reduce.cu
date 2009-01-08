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
using namespace ljgpu::gpu::reduce;

namespace ljgpu { namespace gpu
{

/**
 * blockwise sum
 */
template<typename input_type, typename output_type>
__global__ void sum(input_type const* g_in, output_type* g_block_sum, uint n)
{
    extern __shared__ output_type s_sum[];

    // load values from global device memory
    output_type sum = 0;
    for (uint i = GTID; i < n; i += GTDIM) {
	sum += g_in[i];
    }
    // sum for this thread
    s_sum[TID] = sum;
    __syncthreads();

    // compute sum for all threads in block
    if (TID < 256) {
	sum = sum + s_sum[TID + 256];
	s_sum[TID] = sum;
    }
    __syncthreads();
    if (TID < 128) {
	sum = sum + s_sum[TID + 128];
	s_sum[TID] = sum;
    }
    __syncthreads();
    if (TID < 64) {
	sum = sum + s_sum[TID + 64];
	s_sum[TID] = sum;
    }
    __syncthreads();
    if (TID < 32) {
	sum = sum + s_sum[TID + 32];
	s_sum[TID] = sum;
    }
    // no further syncs needed within execution warp of 32 threads
    if (TID < 16) {
	sum = sum + s_sum[TID + 16];
	s_sum[TID] = sum;
    }
    if (TID < 8) {
	sum = sum + s_sum[TID + 8];
	s_sum[TID] = sum;
    }
    if (TID < 4) {
	sum = sum + s_sum[TID + 4];
	s_sum[TID] = sum;
    }
    if (TID < 2) {
	sum = sum + s_sum[TID + 2];
	s_sum[TID] = sum;
    }
    if (TID < 1) {
	sum = sum + s_sum[TID + 1];
	// store block sum in global memory
	g_block_sum[blockIdx.x] = sum;
    }
}

/**
 * blockwise absolute maximum
 */
template <typename input_type, typename coalesced_input_type, typename output_type>
__global__ void max(coalesced_input_type const* g_in, output_type* g_block_max, uint n)
{
    extern __shared__ output_type s_vv[];

    // load vectors from global device memory
    output_type vv = 0;
    for (uint i = GTID; i < n; i += GTDIM) {
	input_type v = unpack(g_in[i]);
	vv = fmaxf(vv, v * v);
    }
    // maximum squared value for this thread
    s_vv[TID] = vv;
    __syncthreads();

    // compute maximum for all threads in block
    if (TID < 256) {
	vv = fmaxf(vv, s_vv[TID + 256]);
	s_vv[TID] = vv;
    }
    __syncthreads();
    if (TID < 128) {
	vv = fmaxf(vv, s_vv[TID + 128]);
	s_vv[TID] = vv;
    }
    __syncthreads();
    if (TID < 64) {
	vv = fmaxf(vv, s_vv[TID + 64]);
	s_vv[TID] = vv;
    }
    __syncthreads();
    if (TID < 32) {
	vv = fmaxf(vv, s_vv[TID + 32]);
	s_vv[TID] = vv;
    }
    // no further syncs needed within execution warp of 32 threads
    if (TID < 16) {
	vv = fmaxf(vv, s_vv[TID + 16]);
	s_vv[TID] = vv;
    }
    if (TID < 8) {
	vv = fmaxf(vv, s_vv[TID + 8]);
	s_vv[TID] = vv;
    }
    if (TID < 4) {
	vv = fmaxf(vv, s_vv[TID + 4]);
	s_vv[TID] = vv;
    }
    if (TID < 2) {
	vv = fmaxf(vv, s_vv[TID + 2]);
	s_vv[TID] = vv;
    }
    if (TID < 1) {
	vv = fmaxf(vv, s_vv[TID + 1]);
	// store block absolute maximum in global memory
	g_block_max[blockIdx.x] = sqrtf(vv);
    }
}

/**
 * device function wrappers
 */
cuda::function<void(float const*, dfloat*, uint)>
	       reduce::sum(gpu::sum);
cuda::function<void(float4 const*, float*, uint),
	       void(float2 const*, float*, uint)>
	       reduce::max(gpu::max<float3>, gpu::max<float2>);

}} // namespace ljgpu::gpu
