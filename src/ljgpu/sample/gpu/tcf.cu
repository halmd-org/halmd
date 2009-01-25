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
#include <ljgpu/math/gpu/accum.cuh>
#include <ljgpu/math/gpu/vector2d.cuh>
#include <ljgpu/math/gpu/vector3d.cuh>
#include <ljgpu/sample/gpu/tcf.hpp>

namespace ljgpu { namespace cu { namespace tcf
{

enum { THREADS = gpu::tcf_base::THREADS };
enum { WARP_SIZE = gpu::tcf_base::WARP_SIZE };

/** block count */
__constant__ unsigned int* g_n;
/** block mean value */
__constant__ dfloat* g_m;
/** block variance */
__constant__ dfloat* g_v;

template <int threads>
__device__ void reduce(unsigned int& n, dfloat& m, dfloat& v, unsigned int s_n[], dfloat s_m[], dfloat s_v[])
{
    if (TID < threads) {
	accumulator::add(n, m, v, s_n[TID + threads], s_m[TID + threads], s_v[TID + threads]);
	s_n[TID] = n;
	s_m[TID] = m;
	s_v[TID] = v;
    }
    // no further syncs needed within execution warp of 32 threads
    if (threads >= WARP_SIZE) __syncthreads();

    reduce<threads / 2>(n, m, v, s_n, s_m, s_v);
}

template <>
__device__ void reduce<1>(unsigned int& n, dfloat& m, dfloat& v, unsigned int s_n[], dfloat s_m[], dfloat s_v[])
{
    if (TID < 1) {
	accumulator::add(n, m, v, s_n[TID + 1], s_m[TID + 1], s_v[TID + 1]);
    }
}

template <typename vector_type,
	  dfloat (*correlation_function)(vector_type const&, vector_type const&),
	  typename coalesced_vector_type>
__global__ void accumulate(coalesced_vector_type const* g_in, coalesced_vector_type const* g_in0, uint n)
{
    __shared__ unsigned int s_n[THREADS];
    __shared__ dfloat s_m[THREADS];
    __shared__ dfloat s_v[THREADS];

    unsigned int count = 0;
    dfloat mean = 0;
    dfloat variance = 0;

    // load values from global device memory
    for (uint i = GTID; i < n; i += GTDIM) {
	accumulator::add(count, mean, variance, correlation_function(g_in[i], g_in0[i]));
    }
    // reduced value for this thread
    s_n[TID] = count;
    s_m[TID] = mean;
    s_v[TID] = variance;
    __syncthreads();

    // compute reduced value for all threads in block
    reduce<THREADS / 2>(count, mean, variance, s_n, s_m, s_v);

    if (TID < 1) {
	// store block reduced value in global memory
	g_n[blockIdx.x] = count;
	g_m[blockIdx.x] = mean;
	g_v[blockIdx.x] = variance;
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
cuda::function<void (float4 const*, float4 const*, uint)>
    tcf<3>::mean_square_displacement(cu::tcf::accumulate<cu::vector<float, 3>, cu::tcf::mean_square_displacement>);
cuda::function<void (float4 const*, float4 const*, uint)>
    tcf<3>::mean_quartic_displacement(cu::tcf::accumulate<cu::vector<float, 3>, cu::tcf::mean_quartic_displacement>);
cuda::function<void (float4 const*, float4 const*, uint)>
    tcf<3>::velocity_autocorrelation(cu::tcf::accumulate<cu::vector<float, 3>, cu::tcf::velocity_autocorrelation>);

cuda::function<void (float2 const*, float2 const*, uint)>
    tcf<2>::mean_square_displacement(cu::tcf::accumulate<cu::vector<float, 2>, cu::tcf::mean_square_displacement>);
cuda::function<void (float2 const*, float2 const*, uint)>
    tcf<2>::mean_quartic_displacement(cu::tcf::accumulate<cu::vector<float, 2>, cu::tcf::mean_quartic_displacement>);
cuda::function<void (float2 const*, float2 const*, uint)>
    tcf<2>::velocity_autocorrelation(cu::tcf::accumulate<cu::vector<float, 2>, cu::tcf::velocity_autocorrelation>);

/**
 * device symbol wrappers
 */
cuda::symbol<unsigned int*> tcf_base::count(cu::tcf::g_n);
cuda::symbol<dfloat*> tcf_base::mean(cu::tcf::g_m);
cuda::symbol<dfloat*> tcf_base::variance(cu::tcf::g_v);

}} // namespace ljgpu::gpu
