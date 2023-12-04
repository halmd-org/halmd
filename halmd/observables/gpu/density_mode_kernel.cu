/*
 * Copyright © 2008-2019  Felix Höfling
 * Copyright © 2021       Jaslo Ziska
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

#include <halmd/algorithm/gpu/reduction.cuh>
#include <halmd/algorithm/gpu/transform.cuh>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/observables/gpu/density_mode_kernel.hpp>
#include <halmd/utility/gpu/thread.cuh>

using namespace halmd::algorithm::gpu;

namespace halmd {
namespace observables {
namespace gpu {
namespace density_mode_kernel {

// FIXME provide complex data type for CUDA

/**
 *  compute exp(i q·r) for each particle/wavevector pair
 *  and sum results wavevector-wise within a block
 *
 *  @returns block sums of sin(q·r), cos(q·r) for each wavevector
 */
template <int dimension>
__global__ void compute(
    cudaTextureObject_t t_wavevector
  , float4 const* g_r
  , unsigned int const* g_idx, int npart
  , float* g_sin_block, float* g_cos_block, int nq
)
{
    typedef fixed_vector<float, dimension> vector_type;
    typedef typename density_mode_wrapper<dimension>::coalesced_vector_type coalesced_vector_type;

    float sin_;
    float cos_;

    // outer loop over wavevectors
    for (int i=0; i < nq; i++) {
        vector_type q = tex1Dfetch<coalesced_vector_type>(t_wavevector, i);
        sin_ = 0;
        cos_ = 0;
        for (int j = GTID; j < npart; j += GTDIM) {
            // retrieve particle position via index array
            unsigned int idx = g_idx[j];
            vector_type r = g_r[idx];

            float q_r = inner_prod(q, r);
            sin_ += sin(q_r);
            cos_ += cos(q_r);
        }

        // accumulate results within block
        reduce<sum_>(sin_);
        reduce<sum_>(cos_);

        if (TID == 0) {
            g_sin_block[i * BDIM + BID] = sin_;
            g_cos_block[i * BDIM + BID] = cos_;
        }
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
    float s_sum;
    float c_sum;

    // outer loop over wavevectors, distributed over block grid
    for (int i = BID; i < nq; i += BDIM) {
        s_sum = 0;
        c_sum = 0;
        for (int j = TID; j < bdim; j += TDIM) {
            s_sum += g_sin_block[i * bdim + j];
            c_sum += g_cos_block[i * bdim + j];
        }

        // accumulate results within block
        reduce<sum_>(s_sum);
        reduce<sum_>(c_sum);

        // store result in global memory
        if (TID == 0) {
            g_sin[i] = s_sum;
            g_cos[i] = c_sum;
        }
    }
}

} // namespace density_mode_kernel

template <int dimension>
density_mode_wrapper<dimension> density_mode_wrapper<dimension>::kernel = {
    density_mode_kernel::compute<dimension>
  , density_mode_kernel::finalise
};

template class density_mode_wrapper<3>;
template class density_mode_wrapper<2>;

} // namespace gpu
} // namespace observables
} // namespace halmd
