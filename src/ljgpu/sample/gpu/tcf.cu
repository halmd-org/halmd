/* Time correlation functions for CUDA
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
#include <ljgpu/math/gpu/vector2d.cuh>
#include <ljgpu/math/gpu/vector3d.cuh>
#include <ljgpu/sample/gpu/tcf.hpp>
using namespace ljgpu::gpu;

namespace ljgpu { namespace cu { namespace tcf
{

enum { THREADS = tcf_base::THREADS };
enum { WARP_SIZE = tcf_base::WARP_SIZE };

template <int threads>
__device__ void reduce(accumulator<dfloat>& acc, unsigned int s_n[], dfloat s_m[], dfloat s_v[])
{
    if (TID < threads) {
	acc += accumulator<dfloat>(s_n[TID + threads], s_m[TID + threads], s_v[TID + threads]);
	if (threads > 1) {
	    s_m[TID] = acc.mean();
	    s_v[TID] = acc.var();
	    s_n[TID] = acc.count();
	}
    }
    // no further syncs needed within execution warp of 32 threads
    if (threads >= WARP_SIZE) __syncthreads();

    reduce<threads / 2>(acc, s_n, s_m, s_v);
}
template <>
__device__ void reduce<1>(accumulator<dfloat>& acc, unsigned int s_n[], dfloat s_m[], dfloat s_v[])
{
    if (TID < 1) {
	acc += accumulator<dfloat>(s_n[TID + 1], s_m[TID + 1], s_v[TID + 1]);
    }
}

template <typename vector_type,
	  dfloat (*correlation_function)(vector_type const&, vector_type const&),
	  typename coalesced_vector_type>
__global__ void accumulate(coalesced_vector_type const* g_in, coalesced_vector_type const* g_in0, accumulator<dfloat>* g_result, uint n)
{
    __shared__ unsigned int s_n[THREADS];
    __shared__ dfloat s_m[THREADS];
    __shared__ dfloat s_v[THREADS];

    // load values from global device memory
    accumulator<dfloat> acc;
    for (uint i = GTID; i < n; i += GTDIM) {
	acc += correlation_function(g_in[i], g_in0[i]);
    }
    // reduced value for this thread
    s_n[TID] = acc.count();
    s_m[TID] = acc.mean();
    s_v[TID] = acc.var();
    __syncthreads();

    // compute reduced value for all threads in block
    reduce<THREADS / 2>(acc, s_n, s_m, s_v);

    if (TID < 1) {
	// store block reduced value in global memory
	g_result[blockIdx.x] = acc;
    }
}

template <typename vector_type>
__device__ dfloat mean_square_displacement(vector_type const& r, vector_type const& r0)
{
    vector_type const dr = r - r0;
    return dr * dr;
}

template <typename vector_type>
__device__ dfloat mean_quartic_displacement(vector_type const& r, vector_type const& r0)
{
    vector_type const dr = r - r0;
    dfloat const rr = dr * dr;
    return rr * rr;
}

template <typename vector_type>
__device__ dfloat velocity_autocorrelation(vector_type const& v, vector_type const& v0)
{
    return v * v0;
}

}}} // namespace ljgpu::cu::tcf

namespace ljgpu { namespace gpu
{

/**
 * device function wrappers
 */
cuda::function<void (float4 const*, float4 const*, accumulator<dfloat>*, uint)>
    tcf<3>::mean_square_displacement(cu::tcf::accumulate<cu::vector<float, 3>, cu::tcf::mean_square_displacement>);
cuda::function<void (float4 const*, float4 const*, accumulator<dfloat>*, uint)>
    tcf<3>::mean_quartic_displacement(cu::tcf::accumulate<cu::vector<float, 3>, cu::tcf::mean_quartic_displacement>);
cuda::function<void (float4 const*, float4 const*, accumulator<dfloat>*, uint)>
    tcf<3>::velocity_autocorrelation(cu::tcf::accumulate<cu::vector<float, 3>, cu::tcf::velocity_autocorrelation>);

cuda::function<void (float2 const*, float2 const*, accumulator<dfloat>*, uint)>
    tcf<2>::mean_square_displacement(cu::tcf::accumulate<cu::vector<float, 2>, cu::tcf::mean_square_displacement>);
cuda::function<void (float2 const*, float2 const*, accumulator<dfloat>*, uint)>
    tcf<2>::mean_quartic_displacement(cu::tcf::accumulate<cu::vector<float, 2>, cu::tcf::mean_quartic_displacement>);
cuda::function<void (float2 const*, float2 const*, accumulator<dfloat>*, uint)>
    tcf<2>::velocity_autocorrelation(cu::tcf::accumulate<cu::vector<float, 2>, cu::tcf::velocity_autocorrelation>);

}} // namespace ljgpu::gpu
