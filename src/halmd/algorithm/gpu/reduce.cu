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
#include <halmd/algorithm/gpu/reduce.hpp>
#include <halmd/math/gpu/vector2d.cuh>
#include <halmd/math/gpu/vector3d.cuh>

namespace halmd { namespace cu { namespace algorithm
{

enum { THREADS = gpu::reduce::THREADS };

/**
 * parallel reduction
 */
template <typename input_type, typename output_type,
          typename input_transform,
          typename reduce_transform,
          typename output_transform,
          typename coalesced_input_type, typename coalesced_output_type>
__device__ void accumulate(coalesced_input_type const* g_in, coalesced_output_type* g_block_sum, uint n)
{
    __shared__ output_type s_vv[THREADS];

    // load values from global device memory
    output_type vv = 0;
    for (uint i = GTID; i < n; i += GTDIM) {
        output_type v = transform<input_transform, input_type, output_type>(g_in[i]);
        vv = transform<reduce_transform>(vv, v);
    }
    // reduced value for this thread
    s_vv[TID] = vv;
    __syncthreads();

    // compute reduced value for all threads in block
    reduce<THREADS / 2, reduce_transform>(vv, s_vv);

    if (TID < 1) {
        // store block reduced value in global memory
        g_block_sum[blockIdx.x] = transform<output_transform, output_type, output_type>(vv);
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
    accumulate<input_type, output_type, square_, max_, sqrt_>(g_in, g_block_max, n);
}

}}} // namespace halmd::cu::algorithm

namespace halmd { namespace gpu
{

/**
 * device function wrappers
 */
cuda::function<void(float const*, dsfloat*, uint),
               void(float4 const*, float4*, uint),
               void(float2 const*, float2*, uint)>
    reduce::sum(cu::algorithm::sum<float, dsfloat>,
                cu::algorithm::sum<cu::vector<float, 3>, cu::vector<float, 3> >,
                cu::algorithm::sum<cu::vector<float, 2>, cu::vector<float, 2> >);
cuda::function<void(float4 const*, dsfloat*, uint),
               void(float2 const*, dsfloat*, uint)>
    reduce::sum_of_squares(cu::algorithm::sum_of_squares<cu::vector<float, 3>, dsfloat>,
                           cu::algorithm::sum_of_squares<cu::vector<float, 2>, dsfloat>);
cuda::function<void(float4 const*, float*, uint),
               void(float2 const*, float*, uint)>
    reduce::max(cu::algorithm::max<cu::vector<float, 3>, float>,
                cu::algorithm::max<cu::vector<float, 2>, float>);

}} // namespace halmd::gpu
