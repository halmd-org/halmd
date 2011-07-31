/*
 * Copyright © 2008-2010  Peter Colberg
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

#include <boost/preprocessor/repetition/enum_params.hpp>

#include <halmd/algorithm/gpu/reduce_kernel.hpp>
#include <halmd/algorithm/gpu/reduction.cuh>
#include <halmd/utility/gpu/thread.cuh>

namespace halmd {
namespace algorithm {
namespace gpu {
namespace reduce_kernel {

/**
 * parallel reduction
 */
template <
    typename reduce_transform
  , typename input_type
  , typename coalesced_input_type
  , typename output_type
  , typename coalesced_output_type
  , typename input_transform
  , typename output_transform
  , int threads
>
__global__ void reduce(coalesced_input_type const* g_in, coalesced_output_type* g_block_sum, uint n)
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
    gpu::reduce<threads / 2, reduce_transform>(vv, s_vv);

    if (TID < 1) {
        // store block reduced value in global memory
        g_block_sum[blockIdx.x] = transform<output_transform, output_type, output_type>(vv);
    }
}

} // namespace reduce_kernel

//
// To avoid repeating template arguments and at the same time not
// define an ugly macro, we use BOOST_PP_ENUM_PARAMS to generate
// the template arguments. The meaningless names (T0, T1, …) will
// never show up in compile error messages, as the compiler uses
// the template argument names of the *declaration*.
//

template <BOOST_PP_ENUM_PARAMS(7, typename T)>
reduce_wrapper<BOOST_PP_ENUM_PARAMS(7, T)> const reduce_wrapper<BOOST_PP_ENUM_PARAMS(7, T)>::kernel = {
    reduce_kernel::reduce<BOOST_PP_ENUM_PARAMS(7, T), 512>
  , reduce_kernel::reduce<BOOST_PP_ENUM_PARAMS(7, T), 256>
  , reduce_kernel::reduce<BOOST_PP_ENUM_PARAMS(7, T), 128>
  , reduce_kernel::reduce<BOOST_PP_ENUM_PARAMS(7, T),  64>
  , reduce_kernel::reduce<BOOST_PP_ENUM_PARAMS(7, T),  32>
};

} // namespace algorithm
} // namespace gpu
} // namespace halmd
