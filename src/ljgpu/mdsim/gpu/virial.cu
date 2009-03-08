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
#include <ljgpu/mdsim/gpu/virial.cuh>
#include <ljgpu/mdsim/gpu/virial.hpp>

namespace ljgpu { namespace cu { namespace virial
{

enum { THREADS = gpu::virial::THREADS };

/**
 * parallel reduction
 */
__global__ void sum(float4 const* g_virial, float4 const* g_v,
		    dfloat* g_block_sum, uint n)
{
    __shared__ dfloat s_v0[THREADS];
    __shared__ dfloat s_v1[THREADS];
    __shared__ dfloat s_v2[THREADS];
    __shared__ dfloat s_v3[THREADS];

    dfloat v0 = 0;
    dfloat v1 = 0;
    dfloat v2 = 0;
    dfloat v3 = 0;

    // load values from global device memory
    for (uint i = GTID; i < n; i += GTDIM) {
	vector<float, 3> v = g_v[i];
	vector<float, 4> virial = tensor(v * v, v) + g_virial[i];
	v0 += virial.x;
	v1 += virial.y;
	v2 += virial.z;
	v3 += virial.w;
    }
    // reduced value for this thread
    s_v0[TID] = v0;
    s_v1[TID] = v1;
    s_v2[TID] = v2;
    s_v3[TID] = v3;
    __syncthreads();

    // compute reduced value for all threads in block
    reduce<THREADS / 2, quaternion_sum_>(v0, v1, v2, v3, s_v0, s_v1, s_v2, s_v3);

    if (TID < 1) {
	// store block reduced value in global memory
	g_block_sum[blockIdx.x] = v0;
	g_block_sum[blockIdx.x + gridDim.x] = v1;
	g_block_sum[blockIdx.x + 2 * gridDim.x] = v2;
	g_block_sum[blockIdx.x + 3 * gridDim.x] = v3;
    }
}

__global__ void sum(float2 const* g_virial, float2 const* g_v,
		    dfloat* g_block_sum, uint n)
{
    __shared__ dfloat s_v0[THREADS];
    __shared__ dfloat s_v1[THREADS];

    dfloat v0 = 0;
    dfloat v1 = 0;

    // load values from global device memory
    for (uint i = GTID; i < n; i += GTDIM) {
	vector<float, 2> v = g_v[i];
	vector<float, 2> virial = tensor(v * v, v) + g_virial[i];
	v0 += virial.x;
	v1 += virial.y;
    }
    // reduced value for this thread
    s_v0[TID] = v0;
    s_v1[TID] = v1;
    __syncthreads();

    // compute reduced value for all threads in block
    reduce<THREADS / 2, complex_sum_>(v0, v1, s_v0, s_v1);

    if (TID < 1) {
	// store block reduced value in global memory
	g_block_sum[blockIdx.x] = v0;
	g_block_sum[blockIdx.x + gridDim.x] = v1;
    }
}

__global__ void sum(float4 const* g_virial, float4 const* g_v, uint const* g_tag,
		    dfloat* g_block_sum, uint n, uint mpart)
{
    __shared__ dfloat s_v0[THREADS];
    __shared__ dfloat s_v1[THREADS];
    __shared__ dfloat s_v2[THREADS];
    __shared__ dfloat s_v3[THREADS];

    dfloat v0 = 0;
    dfloat v1 = 0;
    dfloat v2 = 0;
    dfloat v3 = 0;
    dfloat v4 = 0;
    dfloat v5 = 0;
    dfloat v6 = 0;
    dfloat v7 = 0;

    // load values from global device memory
    for (uint i = GTID; i < n; i += GTDIM) {
	vector<float, 3> v = g_v[i];
	vector<float, 4> virial = tensor(v * v, v) + g_virial[i];
	if (g_tag[i] < mpart) {
	    v0 += virial.x;
	    v1 += virial.y;
	    v2 += virial.z;
	    v3 += virial.w;
	}
	else {
	    v4 += virial.x;
	    v5 += virial.y;
	    v6 += virial.z;
	    v7 += virial.w;
	}
    }
    // reduced value for this thread
    s_v0[TID] = v0;
    s_v1[TID] = v1;
    s_v2[TID] = v2;
    s_v3[TID] = v3;
    __syncthreads();

    // compute reduced value for all threads in block
    reduce<THREADS / 2, quaternion_sum_>(v0, v1, v2, v3, s_v0, s_v1, s_v2, s_v3);

    if (TID < 1) {
	// store block reduced value in global memory
	g_block_sum[blockIdx.x] = v0;
	g_block_sum[blockIdx.x + gridDim.x] = v1;
	g_block_sum[blockIdx.x + 2 * gridDim.x] = v2;
	g_block_sum[blockIdx.x + 3 * gridDim.x] = v3;
    }

    // reduced value for this thread
    __syncthreads();
    s_v0[TID] = v4;
    s_v1[TID] = v5;
    s_v2[TID] = v6;
    s_v3[TID] = v7;
    __syncthreads();

    // compute reduced value for all threads in block
    reduce<THREADS / 2, quaternion_sum_>(v4, v5, v6, v7, s_v0, s_v1, s_v2, s_v3);

    if (TID < 1) {
	// store block reduced value in global memory
	g_block_sum[blockIdx.x + 4 * gridDim.x] = v4;
	g_block_sum[blockIdx.x + 5 * gridDim.x] = v5;
	g_block_sum[blockIdx.x + 6 * gridDim.x] = v6;
	g_block_sum[blockIdx.x + 7 * gridDim.x] = v7;
    }
}

__global__ void sum(float2 const* g_virial, float2 const* g_v, uint const* g_tag,
		    dfloat* g_block_sum, uint n, uint mpart)
{
    __shared__ dfloat s_v0[THREADS];
    __shared__ dfloat s_v1[THREADS];
    __shared__ dfloat s_v2[THREADS];
    __shared__ dfloat s_v3[THREADS];

    dfloat v0 = 0;
    dfloat v1 = 0;
    dfloat v2 = 0;
    dfloat v3 = 0;

    // load values from global device memory
    for (uint i = GTID; i < n; i += GTDIM) {
	vector<float, 2> v = g_v[i];
	vector<float, 2> virial = tensor(v * v, v) + g_virial[i];
	if (g_tag[i] < mpart) {
	    v0 += virial.x;
	    v1 += virial.y;
	}
	else {
	    v2 += virial.x;
	    v3 += virial.y;
	}
    }
    // reduced value for this thread
    s_v0[TID] = v0;
    s_v1[TID] = v1;
    s_v2[TID] = v2;
    s_v3[TID] = v3;
    __syncthreads();

    // compute reduced value for all threads in block
    reduce<THREADS / 2, quaternion_sum_>(v0, v1, v2, v3, s_v0, s_v1, s_v2, s_v3);

    if (TID < 1) {
	// store block reduced value in global memory
	g_block_sum[blockIdx.x] = v0;
	g_block_sum[blockIdx.x + gridDim.x] = v1;
	g_block_sum[blockIdx.x + 2 * gridDim.x] = v2;
	g_block_sum[blockIdx.x + 3 * gridDim.x] = v3;
    }
}

}}} // namespace ljgpu::cu::virial

namespace ljgpu { namespace gpu
{

/**
 * device function wrappers
 */
cuda::function<
    void(float4 const*, float4 const*, dfloat*, uint),
    void(float2 const*, float2 const*, dfloat*, uint),
    void(float4 const*, float4 const*, uint const*, dfloat*, uint, uint),
    void(float2 const*, float2 const*, uint const*, dfloat*, uint, uint)>
    virial::sum(cu::virial::sum, cu::virial::sum, cu::virial::sum, cu::virial::sum);

}} // namespace ljgpu::gpu
