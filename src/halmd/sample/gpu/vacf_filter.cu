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

#include <halmd/algorithm/gpu/base.cuh>
#include <halmd/algorithm/gpu/reduce.cuh>
#include <halmd/math/gpu/vector2d.cuh>
#include <halmd/math/gpu/vector3d.cuh>
#include <halmd/sample/gpu/vacf_filter.hpp>

namespace halmd { namespace cu { namespace vacf_filter
{

enum { BLOCKS = halmd::gpu::vacf_filter<>::BLOCKS };
enum { THREADS = halmd::gpu::vacf_filter<>::THREADS };

/**
 * calculate block maximum squared displacement, and velocity autocorrelations
 */
template <typename vector_type, typename coalesced_vector_type>
__global__ void accumulate(
    coalesced_vector_type const* g_r
  , coalesced_vector_type const* g_r0
  , coalesced_vector_type const* g_v
  , coalesced_vector_type const* g_v0
  , uint npart
  , float* g_vac //< velocity autocorrelation per particle
  , float* g_msd //< maximum squared displacement per block
)
{
    enum { dimension = vector_type::static_size };
    __shared__ float s_msd[THREADS];
    float msd = 0;

    for (uint i = GTID; i < npart; i += GTDIM) {
        vector_type r = g_r[i], r0 = g_r0[i];
        vector_type v = g_v[i], v0 = g_v0[i];
        g_vac[i] = v * v0;
        vector_type dr = r - r0;
        msd = transform<max_>(msd, dr * dr);
    }

    // reduced value for this thread
    s_msd[TID] = msd;
    __syncthreads();

    // compute reduced value for all threads in block
    reduce<THREADS / 2, max_>(msd, s_msd);

    if (TID < 1) {
        // store block reduced value in global memory
        g_msd[blockIdx.x] = msd;
    }
}

/**
 * assign sort index based on squared displacement relative to maximum
 */
template <typename vector_type, typename coalesced_vector_type>
__global__ void assign_index_from_msd(
    coalesced_vector_type const* g_r
  , coalesced_vector_type const* g_r0
  , uint npart
  , float const* g_msd
  , unsigned int* g_index
)
{
    enum { dimension = vector_type::static_size };
    __shared__ float s_msd[BLOCKS];
    float msd = 0;

    // compute maximum squared displacement from block reduced values
    for (uint i = TID; i < BLOCKS; i += TDIM) {
        s_msd[i] = g_msd[i];
    }
    __syncthreads();
    for (uint i = 0; i < BLOCKS; ++i) {
        msd = transform<max_>(msd, s_msd[i]);
    }

    for (uint i = GTID; i < npart; i += GTDIM) {
        vector_type r = g_r[i], r0 = g_r0[i];
        vector_type dr = r - r0;
        // FIXME unit test to check rescaling to [0, UINT_MAX]
        g_index[i] = float2uint(__saturatef((dr * dr) / msd) * UINT_MAX, cudaRoundNearest);
    }
}

}}} // namespace halmd::cu::vacf_filter

namespace halmd { namespace gpu
{

/**
 * device function wrappers
 */

cuda::function<void (float4 const*, float4 const*, float4 const*, float4 const*, uint, float*, float*)>
    vacf_filter<3>::accumulate(cu::vacf_filter::accumulate<cu::vector<float, 3> >);
cuda::function<void (float4 const*, float4 const*, uint, float const*, unsigned int*)>
    vacf_filter<3>::assign_index_from_msd(cu::vacf_filter::assign_index_from_msd<cu::vector<float, 3> >);

cuda::function<void (float2 const*, float2 const*, float2 const*, float2 const*, uint, float*, float*)>
    vacf_filter<2>::accumulate(cu::vacf_filter::accumulate<cu::vector<float, 2> >);
cuda::function<void (float2 const*, float2 const*, uint, float const*, unsigned int*)>
    vacf_filter<2>::assign_index_from_msd(cu::vacf_filter::assign_index_from_msd<cu::vector<float, 2> >);

}} // namespace halmd::gpu
