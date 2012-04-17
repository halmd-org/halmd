/*
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
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

#include <boost/utility/enable_if.hpp>

#include <halmd/numeric/blas/blas.hpp>
#include <halmd/observables/gpu/density_mode_kernel.hpp>
#include <halmd/utility/gpu/thread.cuh>
#include <halmd/utility/gpu/variant.cuh>

#define MAX_BLOCK_SIZE 512

using namespace boost;
using namespace halmd::utility::gpu; //< variant, map, pair

namespace halmd {
namespace observables {
namespace gpu {
namespace density_mode_kernel {

// pass wavevectors via texture
texture<variant<map<pair<int_<3>, float4>, pair<int_<2>, float2> > > > q_;

// global constants
__constant__ uint nq_;        // number of wavevectors

// recursive reduction function,
// terminate for threads=0
template <unsigned threads, typename T>
__device__ typename disable_if_c<threads>::type
sum_reduce(T*, T*) {}

// reduce two array simultaneously by summation,
// size of a,b must be at least 2 * threads
template <unsigned threads, typename T>
__device__ typename enable_if_c<threads>::type
sum_reduce(T* a, T* b)
{
    if (TID < threads) {
        a[TID] += a[TID + threads];
        b[TID] += b[TID + threads];
    }
    if (threads >= warpSize) {
        __syncthreads();
    }

    // recursion ends by calling sum_reduce<0>
    sum_reduce<threads / 2>(a, b);
}


/* FIXME
typedef void (*sum_reduce_type)(float*, float*);
__device__ sum_reduce_type sum_reduce_select[] = {
    &sum_reduce<0>, &sum_reduce<1>, &sum_reduce<2>, &sum_reduce<4>,
    &sum_reduce<8>, &sum_reduce<16>, &sum_reduce<32>, &sum_reduce<64>,
    &sum_reduce<128>, &sum_reduce<256>
};
*/

// FIXME provide complex data type for CUDA

/**
 *  compute exp(i q·r) for each particle/wavevector pair
 *  and sum results wavevector-wise within a block
 *
 *  @returns block sums of sin(q·r), cos(q·r) for each wavevector
 */
template <typename vector_type, typename coalesced_vector_type>
__global__ void compute(coalesced_vector_type const* g_r, uint npart, float* g_sin_block, float* g_cos_block)
{
    enum { dimension = vector_type::static_size };

    __shared__ float sin_[MAX_BLOCK_SIZE];
    __shared__ float cos_[MAX_BLOCK_SIZE];

    // outer loop over wavevectors
    for (uint i=0; i < nq_; i++) {
        vector_type q = tex1Dfetch(get<dimension>(q_), i);
        sin_[TID] = 0;
        cos_[TID] = 0;
        for (uint j = GTID; j < npart; j += GTDIM) {
            // retrieve particle position
            vector_type r = g_r[j];

            float q_r = inner_prod(q, r);
            sin_[TID] += sin(q_r);
            cos_[TID] += cos(q_r);
        }
        __syncthreads();

        // accumulate results within block
        if (TDIM == 512) sum_reduce<256>(sin_, cos_);
        else if (TDIM == 256) sum_reduce<128>(sin_, cos_);
        else if (TDIM == 128) sum_reduce<64>(sin_, cos_);
        else if (TDIM == 64) sum_reduce<32>(sin_, cos_);
        else if (TDIM == 32) sum_reduce<16>(sin_, cos_);
        else if (TDIM == 16) sum_reduce<8>(sin_, cos_);
        else if (TDIM == 8) sum_reduce<4>(sin_, cos_);

        if (TID == 0) {
            g_sin_block[i * BDIM + BID] = sin_[0];
            g_cos_block[i * BDIM + BID] = cos_[0];
        }
        __syncthreads();
    }
}

/**
 *  reduce block sums for each wavevector separately
 *
 *  @param bdim  number of blocks (grid size) in the preceding call to compute()
 */
__global__ void finalise(
    float const* g_sin_block, float const* g_cos_block
  , float* g_sin, float* g_cos
  , uint bdim)
{
    __shared__ float s_sum[MAX_BLOCK_SIZE];
    __shared__ float c_sum[MAX_BLOCK_SIZE];

    // outer loop over wavevectors, distributed over block grid
    for (uint i = BID; i < nq_; i += BDIM) {
        s_sum[TID] = 0;
        c_sum[TID] = 0;
        for (uint j = TID; j < bdim; j += TDIM) {
            s_sum[TID] += g_sin_block[i * bdim + j];
            c_sum[TID] += g_cos_block[i * bdim + j];
        }
        __syncthreads();

        // accumulate results within block
        if (TDIM == 512) sum_reduce<256>(s_sum, c_sum);
        else if (TDIM == 256) sum_reduce<128>(s_sum, c_sum);
        else if (TDIM == 128) sum_reduce<64>(s_sum, c_sum);
        else if (TDIM == 64) sum_reduce<32>(s_sum, c_sum);
        else if (TDIM == 32) sum_reduce<16>(s_sum, c_sum);
        else if (TDIM == 16) sum_reduce<8>(s_sum, c_sum);
        else if (TDIM == 8) sum_reduce<4>(s_sum, c_sum);

        // store result in global memory
        if (TID == 0) {
            g_sin[i] = s_sum[0];
            g_cos[i] = c_sum[0];
        }
    }
}

} // namespace density_mode_kernel

template <int dimension>
density_mode_wrapper<dimension> const density_mode_wrapper<dimension>::kernel = {
    get<dimension>(density_mode_kernel::q_)
  , density_mode_kernel::nq_
  , density_mode_kernel::compute<fixed_vector<float, dimension> >
  , density_mode_kernel::finalise
};

template class density_mode_wrapper<3>;
template class density_mode_wrapper<2>;

} // namespace observables
} // namespace gpu
} // namespace halmd
