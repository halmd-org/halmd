/*
 * Copyright © 2008-2019  Felix Höfling
 * Copyright © 2015       Nicolas Höft
 * Copyright © 2008-2011  Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <boost/utility/enable_if.hpp>

#include <halmd/numeric/blas/blas.hpp>
#include <halmd/observables/gpu/density_mode_kernel.hpp>
#include <halmd/utility/gpu/thread.cuh>

#define MAX_BLOCK_SIZE 1024

namespace halmd {
namespace observables {
namespace gpu {
namespace density_mode_kernel {

// pass wavevectors via texture
template<int dimension>
struct wavevector
{
    // instantiate a separate texture for each aligned vector type
    typedef texture<typename density_mode_wrapper<dimension>::coalesced_vector_type> type;
    static type tex_;
};
// instantiate static members
template<int dimension> wavevector<dimension>::type wavevector<dimension>::tex_;

// recursive reduction function,
// terminate for threads=0
template <unsigned threads, typename T>
__device__ typename boost::disable_if_c<threads>::type
sum_reduce(T*, T*) {}

// reduce two array simultaneously by summation,
// size of a,b must be at least 2 * threads
template <unsigned threads, typename T>
__device__ typename boost::enable_if_c<threads>::type
sum_reduce(T* a, T* b)
{
    if (TID < threads) {
        a[TID] += a[TID + threads];
        b[TID] += b[TID + threads];
    }

    if (threads >= warpSize) {
        __syncthreads();
    }
    else {
        // on hardware of compute capability ≥ 7.0 (Volta),
        // warps are no longer guaranteed to be executed in lock-step
#if CUDART_VERSION >= 9000
        // select warp lanes with TID < threads,
        // fix compilation of operator<< for large values of 'threads'
        unsigned mask = (1U << (threads & (warpSize - 1))) - 1;
        __syncwarp(mask);
#else
        __syncthreads();    // only needed if the _hardware_ is Volta or later
#endif
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
__global__ void compute(
    coalesced_vector_type const* g_r
  , unsigned int const* g_idx, int npart
  , float* g_sin_block, float* g_cos_block, int nq
)
{
    enum { dimension = vector_type::static_size };

    __shared__ float sin_[MAX_BLOCK_SIZE];
    __shared__ float cos_[MAX_BLOCK_SIZE];

    // outer loop over wavevectors
    for (int i=0; i < nq; i++) {
        vector_type q = tex1Dfetch(wavevector<dimension>::tex_, i);
        sin_[TID] = 0;
        cos_[TID] = 0;
        for (int j = GTID; j < npart; j += GTDIM) {
            // retrieve particle position via index array
            unsigned int idx = g_idx[j];
            vector_type r = g_r[idx];

            float q_r = inner_prod(q, r);
            sin_[TID] += sin(q_r);
            cos_[TID] += cos(q_r);
        }
        __syncthreads();

        // accumulate results within block
        if (TDIM == 1024) sum_reduce<512>(sin_, cos_);
        else if (TDIM == 512) sum_reduce<256>(sin_, cos_);
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
        __syncthreads();    // FIXME needed here? would __syncwarp() be sufficient?
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
  , int nq, int bdim)
{
    __shared__ float s_sum[MAX_BLOCK_SIZE];
    __shared__ float c_sum[MAX_BLOCK_SIZE];

    // outer loop over wavevectors, distributed over block grid
    for (int i = BID; i < nq; i += BDIM) {
        s_sum[TID] = 0;
        c_sum[TID] = 0;
        for (int j = TID; j < bdim; j += TDIM) {
            s_sum[TID] += g_sin_block[i * bdim + j];
            c_sum[TID] += g_cos_block[i * bdim + j];
        }
        __syncthreads();

        // accumulate results within block
        if (TDIM == 1024) sum_reduce<512>(s_sum, c_sum);
        else if (TDIM == 512) sum_reduce<256>(s_sum, c_sum);
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
        __syncthreads();    // FIXME needed here?
    }
}

} // namespace density_mode_kernel

template <int dimension>
density_mode_wrapper<dimension> const density_mode_wrapper<dimension>::kernel = {
    density_mode_kernel::wavevector<dimension>::tex_
  , density_mode_kernel::compute<fixed_vector<float, dimension> >
  , density_mode_kernel::finalise
};

template class density_mode_wrapper<3>;
template class density_mode_wrapper<2>;

} // namespace gpu
} // namespace observables
} // namespace halmd
