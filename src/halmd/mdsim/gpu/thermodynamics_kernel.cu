/*
 * Copyright © 2010  Felix Höfling and Peter Colberg
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
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>

using namespace halmd::algorithm::gpu;

namespace halmd
{
namespace mdsim { namespace gpu
{
namespace thermodynamics_kernel
{

/**
 * parallel reduction, separately for each particle type
 * as provided by g_r
 */
template <
    typename reduce_transform
  , typename input_type
  , typename coalesced_input_type
  , typename output_type
  , typename coalesced_output_type
  , typename coalesced_vector_type
  , typename input_transform
  , typename output_transform
  , int threads
  , int ntypes
>
__global__ void reduce_types(coalesced_input_type const* g_in,
                             coalesced_output_type* g_block_sum,
                             coalesced_vector_type* g_r, uint npart)
{
    typedef fixed_vector<output_type, ntypes> output_vector_type;
    __shared__ output_vector_type s_vv[threads];

    // load values from global device memory
    output_vector_type vv;
    for (uint k = 0; k < ntypes; k++) {
        vv[k] = 0;
    }
    for (uint i = GTID; i < n; i += GTDIM) {
        // load and transform input value
        output_type v = transform<input_transform, input_type, output_type>(g_in[i]);
        // determine particle type
        vector_type r;
        unsigned int type;
        tie(r, type) = untagged<vector_type>(g_r[GTID]);
        // assert(type < ntypes);
        // merge (add, multiply, ...) with accumulator
        vv[type] = transform<reduce_transform>(vv[type], v);
    }
    // reduced value for this thread
    s_vv[TID] = vv;
    __syncthreads();

    // compute reduced value for all threads in block
    gpu::reduce<threads / 2, reduce_transform>(vv, s_vv);

    if (TID < 1) {
        // store block reduced value in global memory
        // in groups of particle types
        for (uint k = 0; k < ntypes; k++) {
            g_block_sum[BID + k * BDIM] =
                transform<output_transform, output_type, output_type>(vv[k]);
        }
    }
}

} // namespace thermodynamics_kernel

}} // namespace mdsim::gpu

template class reduce_wrapper<
    sum_                        // reduce_transform
  , fixed_vector<float, 3>      // input_type
  , float4                      // coalesced_input_type
  , dsfloat                     // output_type
  , dsfloat                     // coalesced_output_type
  , square_                     // input_transform
>;

template class reduce_wrapper<
    sum_                        // reduce_transform
  , fixed_vector<float, 2>      // input_type
  , float4                      // coalesced_input_type
  , dsfloat                     // output_type
  , dsfloat                     // coalesced_output_type
  , square_                     // input_transform
>;

template class reduce_wrapper<
    sum_                        // reduce_transform
  , fixed_vector<float, 3>      // input_type
  , float4                      // coalesced_input_type
  , fixed_vector<dsfloat, 3>    // output_type
  , fixed_vector<dsfloat, 3>    // coalesced_output_type
>;

template class reduce_wrapper<
    sum_                        // reduce_transform
  , fixed_vector<float, 2>      // input_type
  , float4                      // coalesced_input_type
  , fixed_vector<dsfloat, 2>    // output_type
  , fixed_vector<dsfloat, 2>    // coalesced_output_type
>;

template class reduce_wrapper<
    sum_                        // reduce_transform
  , float                       // input_type
  , float                       // coalesced_input_type
  , dsfloat                     // output_type
  , dsfloat                     // coalesced_output_type
>;

} // namespace halmd
