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
__global__ void accumulate(coalesced_vector_type const* g_in, coalesced_vector_type const* g_in0, uint* g_n, dfloat* g_m, dfloat* g_v, uint n)
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

template <int threads>
__device__ void reduce(dfloat& sum, dfloat s_sum[])
{
    if (TID < threads) {
	sum += s_sum[TID + threads];
	s_sum[TID] = sum;
    }
    // no further syncs needed within execution warp of 32 threads
    if (threads >= WARP_SIZE) __syncthreads();

    reduce<threads / 2>(sum, s_sum);
}

template <>
__device__ void reduce<1>(dfloat& sum, dfloat s_sum[])
{
    if (TID < 1) {
	sum += s_sum[TID + 1];
    }
}

template <typename vector_type,
	  dfloat (*correlation_function)(vector_type const&, vector_type const&, vector_type const&),
	  typename coalesced_vector_type,
	  typename uncoalesced_vector_type>
__global__ void accumulate(coalesced_vector_type const* g_in, coalesced_vector_type const* g_in0, uncoalesced_vector_type const q_vector, dfloat* g_sum, uint n)
{
    __shared__ dfloat s_sum[THREADS];

    dfloat sum = 0;

    // load values from global device memory
    for (uint i = GTID; i < n; i += GTDIM) {
	sum += correlation_function(g_in[i], g_in0[i], q_vector);
    }
    // reduced value for this thread
    s_sum[TID] = sum;
    __syncthreads();

    // compute reduced value for all threads in block
    reduce<THREADS / 2>(sum, s_sum);

    if (TID < 1) {
	// store block reduced value in global memory
	g_sum[blockIdx.x] = sum;
    }
}

template <typename vector_type>
__device__ dfloat incoherent_scattering_function(vector_type const& r, vector_type const& r0, vector_type const& q)
{
    return cosf((r - r0) * q);
}

template <int threads>
__device__ void reduce(vector<dfloat, 2>& sum, dfloat s_imag[], dfloat s_real[])
{
    if (TID < threads) {
	sum += vector<dfloat, 2>(s_real[TID + threads], s_imag[TID + threads]);
	s_real[TID] = sum.x;
	s_imag[TID] = sum.y;
    }
    // no further syncs needed within execution warp of 32 threads
    if (threads >= WARP_SIZE) __syncthreads();

    reduce<threads / 2>(sum, s_real, s_imag);
}

template <>
__device__ void reduce<1>(vector<dfloat, 2>& sum, dfloat s_imag[], dfloat s_real[])
{
    if (TID < 1) {
	sum += vector<dfloat, 2>(s_real[TID + 1], s_imag[TID + 1]);
    }
}

template <typename vector_type,
	  vector<dfloat, 2> (*correlation_function)(vector_type const&, vector_type const&),
	  typename coalesced_vector_type,
	  typename uncoalesced_vector_type>
__global__ void accumulate(coalesced_vector_type const* g_in, uncoalesced_vector_type const q_vector, dfloat* g_real, dfloat* g_imag, uint n)
{
    __shared__ dfloat s_real[THREADS];
    __shared__ dfloat s_imag[THREADS];

    vector<dfloat, 2> sum = 0;

    // load values from global device memory
    for (uint i = GTID; i < n; i += GTDIM) {
	sum += correlation_function(g_in[i], q_vector);
    }
    // reduced value for this thread
    s_real[TID] = sum.x;
    s_imag[TID] = sum.y;
    __syncthreads();

    // compute reduced value for all threads in block
    reduce<THREADS / 2>(sum, s_real, s_imag);

    if (TID < 1) {
	// store block reduced value in global memory
	g_real[blockIdx.x] = sum.x;
	g_imag[blockIdx.x] = sum.y;
    }
}

template <typename vector_type>
__device__ vector<dfloat, 2> coherent_scattering_function(vector_type const& r, vector_type const& q)
{
    dfloat const r_q = r * q;
    return vector<dfloat, 2> (cosf(r_q), sinf(r_q));
}

}}} // namespace ljgpu::cu::tcf

namespace ljgpu { namespace gpu
{

/**
 * device function wrappers
 */
cuda::function<void (float4 const*, float4 const*, uint*, dfloat*, dfloat*, uint)>
    tcf<3>::mean_square_displacement(cu::tcf::accumulate<cu::vector<float, 3>, cu::tcf::mean_square_displacement>);
cuda::function<void (float4 const*, float4 const*, uint*, dfloat*, dfloat*, uint)>
    tcf<3>::mean_quartic_displacement(cu::tcf::accumulate<cu::vector<float, 3>, cu::tcf::mean_quartic_displacement>);
cuda::function<void (float4 const*, float4 const*, uint*, dfloat*, dfloat*, uint)>
    tcf<3>::velocity_autocorrelation(cu::tcf::accumulate<cu::vector<float, 3>, cu::tcf::velocity_autocorrelation>);
cuda::function<void (float4 const*, float4 const*, float3 const, dfloat*, uint)>
    tcf<3>::incoherent_scattering_function(cu::tcf::accumulate<cu::vector<float, 3>, cu::tcf::incoherent_scattering_function>);
cuda::function<void (float4 const*, float3 const, dfloat*, dfloat*, uint)>
    tcf<3>::coherent_scattering_function(cu::tcf::accumulate<cu::vector<float, 3>, cu::tcf::coherent_scattering_function>);

cuda::function<void (float2 const*, float2 const*, uint*, dfloat*, dfloat*, uint)>
    tcf<2>::mean_square_displacement(cu::tcf::accumulate<cu::vector<float, 2>, cu::tcf::mean_square_displacement>);
cuda::function<void (float2 const*, float2 const*, uint*, dfloat*, dfloat*, uint)>
    tcf<2>::mean_quartic_displacement(cu::tcf::accumulate<cu::vector<float, 2>, cu::tcf::mean_quartic_displacement>);
cuda::function<void (float2 const*, float2 const*, uint*, dfloat*, dfloat*, uint)>
    tcf<2>::velocity_autocorrelation(cu::tcf::accumulate<cu::vector<float, 2>, cu::tcf::velocity_autocorrelation>);
cuda::function<void (float2 const*, float2 const*, float2 const, dfloat*, uint)>
    tcf<2>::incoherent_scattering_function(cu::tcf::accumulate<cu::vector<float, 2>, cu::tcf::incoherent_scattering_function>);
cuda::function<void (float2 const*, float2 const, dfloat*, dfloat*, uint)>
    tcf<2>::coherent_scattering_function(cu::tcf::accumulate<cu::vector<float, 2>, cu::tcf::coherent_scattering_function>);

}} // namespace ljgpu::gpu
