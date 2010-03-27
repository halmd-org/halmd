/* Parallel reduction kernel
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
#include <halmd/algorithm/gpu/reduce.cuh>
#include <halmd/math/gpu/dsvector.cuh>
#include <halmd/mdsim/backend/gpu/virial.cuh>
#include <halmd/mdsim/backend/gpu/virial.hpp>

namespace halmd { namespace cu { namespace virial
{

enum { THREADS = gpu::virial::THREADS };

/**
 * Virial stress tensor for three-dimensional monodisperse system
 */
__global__ void sum(float4 const* g_virial, float4 const* g_v,
                    vector<dsfloat, 4>* g_block_sum, uint n)
{
    __shared__ vector<dsfloat, 4> s_virial[THREADS];
    vector<dsfloat, 4> virial = 0;

    // load values from global device memory
    for (uint i = GTID; i < n; i += GTDIM) {
        vector<float, 3> v = g_v[i];
        virial += tensor(v * v, v) + g_virial[i];
    }
    // reduced value for this thread
    s_virial[TID] = virial;
    __syncthreads();

    // compute reduced value for all threads in block
    reduce<THREADS / 2, sum_>(virial, s_virial);

    if (TID < 1) {
        // store block reduced value in global memory
        g_block_sum[blockIdx.x] = virial;
    }
}

/**
 * Virial stress tensor for two-dimensional monodisperse system
 */
__global__ void sum(float2 const* g_virial, float2 const* g_v,
                    vector<dsfloat, 2>* g_block_sum, uint n)
{
    __shared__ vector<dsfloat, 2> s_virial[THREADS];
    vector<dsfloat, 2> virial = 0;

    // load values from global device memory
    for (uint i = GTID; i < n; i += GTDIM) {
        vector<float, 2> v = g_v[i];
        virial += tensor(v * v, v) + g_virial[i];
    }
    // reduced value for this thread
    s_virial[TID] = virial;
    __syncthreads();

    // compute reduced value for all threads in block
    reduce<THREADS / 2, sum_>(virial, s_virial);

    if (TID < 1) {
        // store block reduced value in global memory
        g_block_sum[blockIdx.x] = virial;
    }
}

/**
 * Virial stress tensor for three-dimensional bidisperse system
 */
__global__ void sum(float4 const* g_virial, float4 const* g_v, uint const* g_tag,
                    vector<dsfloat, 4>* g_block_sum, uint n, uint mpart)
{
    __shared__ vector<dsfloat, 4> s_virial[THREADS];
    vector<dsfloat, 4> virial_a = 0;
    vector<dsfloat, 4> virial_b = 0;

    // load values from global device memory
    for (uint i = GTID; i < n; i += GTDIM) {
        vector<float, 3> v = g_v[i];
        vector<float, 4> virial = tensor(v * v, v) + g_virial[i];
        if (g_tag[i] < mpart) {
            virial_a += virial;
        }
        else {
            virial_b += virial;
        }
    }
    // reduced value for this thread
    s_virial[TID] = virial_a;
    __syncthreads();

    // compute reduced value for all threads in block
    reduce<THREADS / 2, sum_>(virial_a, s_virial);

    if (TID < 1) {
        // store block reduced value in global memory
        g_block_sum[blockIdx.x] = virial_a;
    }

    // reduced value for this thread
    __syncthreads();
    s_virial[TID] = virial_b;
    __syncthreads();

    // compute reduced value for all threads in block
    reduce<THREADS / 2, sum_>(virial_b, s_virial);

    if (TID < 1) {
        // store block reduced value in global memory
        g_block_sum[blockIdx.x + gridDim.x] = virial_b;
    }
}

/**
 * Virial stress tensor for two-dimensional bidisperse system
 */
__global__ void sum(float2 const* g_virial, float2 const* g_v, uint const* g_tag,
                    vector<dsfloat, 2>* g_block_sum, uint n, uint mpart)
{
    __shared__ vector<dsfloat, 2> s_virial_a[THREADS];
    __shared__ vector<dsfloat, 2> s_virial_b[THREADS];
    vector<dsfloat, 2> virial_a = 0;
    vector<dsfloat, 2> virial_b = 0;

    // load values from global device memory
    for (uint i = GTID; i < n; i += GTDIM) {
        vector<float, 2> v = g_v[i];
        vector<float, 2> virial = tensor(v * v, v) + g_virial[i];
        if (g_tag[i] < mpart) {
            virial_a += virial;
        }
        else {
            virial_b += virial;
        }
    }
    // reduced value for this thread
    s_virial_a[TID] = virial_a;
    s_virial_b[TID] = virial_b;
    __syncthreads();

    // compute reduced value for all threads in block
    reduce<THREADS / 2, complex_sum_>(virial_a, virial_b, s_virial_a, s_virial_b);

    if (TID < 1) {
        // store block reduced value in global memory
        g_block_sum[blockIdx.x] = virial_a;
        g_block_sum[blockIdx.x + gridDim.x] = virial_b;
    }
}

}}} // namespace halmd::cu::virial

namespace halmd { namespace gpu
{

/**
 * device function wrappers
 */
cuda::function<
    void(float4 const*, float4 const*, cu::vector<dsfloat, 4>*, uint),
    void(float2 const*, float2 const*, cu::vector<dsfloat, 2>*, uint),
    void(float4 const*, float4 const*, uint const*, cu::vector<dsfloat, 4>*, uint, uint),
    void(float2 const*, float2 const*, uint const*, cu::vector<dsfloat, 2>*, uint, uint)>
    virial::sum(cu::virial::sum, cu::virial::sum, cu::virial::sum, cu::virial::sum);

}} // namespace halmd::gpu
