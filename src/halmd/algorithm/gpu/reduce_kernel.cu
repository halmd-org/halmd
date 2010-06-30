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

#include <halmd/algorithm/gpu/reduce_kernel.cuh>
#include <halmd/numeric/gpu/blas/dsfloat.cuh>
#include <halmd/numeric/gpu/blas/vector.cuh>
#include <halmd/utility/gpu/thread.cuh>

using namespace halmd::numeric::gpu::blas;

namespace halmd
{
namespace algorithm { namespace gpu
{

/**
 * parallel reduction
 */
template <
    int threads
  , typename input_type
  , typename output_type
  , typename input_transform
  , typename reduce_transform
  , typename output_transform
  , typename coalesced_input_type
  , typename coalesced_output_type
>
__device__ void accumulate(coalesced_input_type const* g_in, coalesced_output_type* g_block_sum, uint n)
{
    __shared__ output_type s_vv[threads];

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
    reduce<threads / 2, reduce_transform>(vv, s_vv);

    if (TID < 1) {
        // store block reduced value in global memory
        g_block_sum[blockIdx.x] = transform<output_transform, output_type, output_type>(vv);
    }
}

/**
 * blockwise sum
 */
template <
    int threads
  , typename input_type
  , typename output_type
  , typename coalesced_input_type
  , typename coalesced_output_type
>
__global__ void sum(coalesced_input_type const* g_in, coalesced_output_type* g_block_sum, uint n)
{
    accumulate<threads, input_type, output_type, identity_, sum_, identity_>(g_in, g_block_sum, n);
}

/**
 * blockwise sum of squares
 */
template <
    int threads
  , typename input_type
  , typename output_type
  , typename coalesced_input_type
  , typename coalesced_output_type
>
__global__ void sum_of_squares(coalesced_input_type const* g_in, coalesced_output_type* g_block_sum, uint n)
{
    accumulate<threads, input_type, output_type, square_, sum_, identity_>(g_in, g_block_sum, n);
}

/**
 * blockwise absolute maximum
 */
template <
    int threads
  , typename input_type
  , typename output_type
  , typename coalesced_input_type
  , typename coalesced_output_type
>
__global__ void max(coalesced_input_type const* g_in, coalesced_output_type* g_block_max, uint n)
{
    accumulate<threads, input_type, output_type, square_, max_, sqrt_>(g_in, g_block_max, n);
}

/**
 * device function wrappers
 */
template <> cuda::function<void(float4 const*, float4*, uint)>
    reduce_wrapper<512, float4, float4>::sum(gpu::sum<512, vector<float, 3>, vector<float, 3> >);
template <> cuda::function<void(float2 const*, float2*, uint)>
    reduce_wrapper<512, float2, float2>::sum(gpu::sum<512, vector<float, 2>, vector<float, 2> >);
template <> cuda::function<void(float4 const*, float4*, uint)>
    reduce_wrapper<256, float4, float4>::sum(gpu::sum<256, vector<float, 3>, vector<float, 3> >);
template <> cuda::function<void(float2 const*, float2*, uint)>
    reduce_wrapper<256, float2, float2>::sum(gpu::sum<256, vector<float, 2>, vector<float, 2> >);
template <> cuda::function<void(float4 const*, float4*, uint)>
    reduce_wrapper<128, float4, float4>::sum(gpu::sum<128, vector<float, 3>, vector<float, 3> >);
template <> cuda::function<void(float2 const*, float2*, uint)>
    reduce_wrapper<128, float2, float2>::sum(gpu::sum<128, vector<float, 2>, vector<float, 2> >);
template <> cuda::function<void(float4 const*, float4*, uint)>
    reduce_wrapper<64, float4, float4>::sum(gpu::sum<64, vector<float, 3>, vector<float, 3> >);
template <> cuda::function<void(float2 const*, float2*, uint)>
    reduce_wrapper<64, float2, float2>::sum(gpu::sum<64, vector<float, 2>, vector<float, 2> >);

template <> cuda::function<void(float4 const*, dsfloat*, uint)>
    reduce_wrapper<512, float4, dsfloat>::sum_of_squares(gpu::sum_of_squares<512, vector<float, 3>, dsfloat>);
template <> cuda::function<void(float2 const*, dsfloat*, uint)>
    reduce_wrapper<512, float2, dsfloat>::sum_of_squares(gpu::sum_of_squares<512, vector<float, 2>, dsfloat>);
template <> cuda::function<void(float4 const*, dsfloat*, uint)>
    reduce_wrapper<256, float4, dsfloat>::sum_of_squares(gpu::sum_of_squares<256, vector<float, 3>, dsfloat>);
template <> cuda::function<void(float2 const*, dsfloat*, uint)>
    reduce_wrapper<256, float2, dsfloat>::sum_of_squares(gpu::sum_of_squares<256, vector<float, 2>, dsfloat>);
template <> cuda::function<void(float4 const*, dsfloat*, uint)>
    reduce_wrapper<128, float4, dsfloat>::sum_of_squares(gpu::sum_of_squares<128, vector<float, 3>, dsfloat>);
template <> cuda::function<void(float2 const*, dsfloat*, uint)>
    reduce_wrapper<128, float2, dsfloat>::sum_of_squares(gpu::sum_of_squares<128, vector<float, 2>, dsfloat>);
template <> cuda::function<void(float4 const*, dsfloat*, uint)>
    reduce_wrapper<64, float4, dsfloat>::sum_of_squares(gpu::sum_of_squares<64, vector<float, 3>, dsfloat>);
template <> cuda::function<void(float2 const*, dsfloat*, uint)>
    reduce_wrapper<64, float2, dsfloat>::sum_of_squares(gpu::sum_of_squares<64, vector<float, 2>, dsfloat>);

template <> cuda::function<void(float4 const*, float*, uint)>
    reduce_wrapper<512, float4, float>::max(gpu::max<512, vector<float, 3>, float>);
template <> cuda::function<void(float2 const*, float*, uint)>
    reduce_wrapper<512, float2, float>::max(gpu::max<512, vector<float, 2>, float>);
template <> cuda::function<void(float4 const*, float*, uint)>
    reduce_wrapper<256, float4, float>::max(gpu::max<256, vector<float, 3>, float>);
template <> cuda::function<void(float2 const*, float*, uint)>
    reduce_wrapper<256, float2, float>::max(gpu::max<256, vector<float, 2>, float>);
template <> cuda::function<void(float4 const*, float*, uint)>
    reduce_wrapper<128, float4, float>::max(gpu::max<128, vector<float, 3>, float>);
template <> cuda::function<void(float2 const*, float*, uint)>
    reduce_wrapper<128, float2, float>::max(gpu::max<128, vector<float, 2>, float>);
template <> cuda::function<void(float4 const*, float*, uint)>
    reduce_wrapper<64, float4, float>::max(gpu::max<64, vector<float, 3>, float>);
template <> cuda::function<void(float2 const*, float*, uint)>
    reduce_wrapper<64, float2, float>::max(gpu::max<64, vector<float, 2>, float>);

}} // namespace algorithm::gpu

} // namespace halmd
