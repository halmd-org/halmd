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
#include <ljgpu/algorithm/gpu/reduce.cuh>
#include <ljgpu/math/gpu/vector2d.cuh>
#include <ljgpu/math/gpu/vector3d.cuh>
#include <ljgpu/sample/gpu/tcf.hpp>

namespace ljgpu { namespace cu { namespace tcf
{

enum { THREADS = gpu::tcf_base::THREADS };

template <typename vector_type,
          dsfloat (*correlation_function)(vector_type const&, vector_type const&),
          typename coalesced_vector_type>
__global__ void accumulate(coalesced_vector_type const* g_in, coalesced_vector_type const* g_in0, uint* g_n, dsfloat* g_m, dsfloat* g_v, uint n)
{
    __shared__ unsigned int s_n[THREADS];
    __shared__ dsfloat s_m[THREADS];
    __shared__ dsfloat s_v[THREADS];

    unsigned int count = 0;
    dsfloat mean = 0;
    dsfloat variance = 0;

    // load values from global device memory
    for (uint i = GTID; i < n; i += GTDIM) {
        transform<accumulate_>(count, mean, variance, correlation_function(g_in[i], g_in0[i]));
    }
    // reduced value for this thread
    s_n[TID] = count;
    s_m[TID] = mean;
    s_v[TID] = variance;
    __syncthreads();

    // compute reduced value for all threads in block
    reduce<THREADS / 2, accumulate_>(count, mean, variance, s_n, s_m, s_v);

    if (TID < 1) {
        // store block reduced value in global memory
        g_n[blockIdx.x] = count;
        g_m[blockIdx.x] = mean;
        g_v[blockIdx.x] = variance;
    }
}

template <typename vector_type>
__device__ dsfloat mean_square_displacement(vector_type const& r, vector_type const& r0)
{
    vector_type const dr = r - r0;
    return dr * dr;
}

template <typename vector_type>
__device__ dsfloat mean_quartic_displacement(vector_type const& r, vector_type const& r0)
{
    vector_type const dr = r - r0;
    dsfloat const rr = dr * dr;
    return rr * rr;
}

template <typename vector_type>
__device__ dsfloat velocity_autocorrelation(vector_type const& v, vector_type const& v0)
{
    return v * v0;
}

template <typename vector_type,
          dsfloat (*correlation_function)(vector_type const&, vector_type const&, vector_type const&),
          typename coalesced_vector_type,
          typename uncoalesced_vector_type>
__global__ void accumulate(coalesced_vector_type const* g_in, coalesced_vector_type const* g_in0, uncoalesced_vector_type const q_vector, dsfloat* g_sum, uint n)
{
    __shared__ dsfloat s_sum[THREADS];

    dsfloat sum = 0;

    // load values from global device memory
    for (uint i = GTID; i < n; i += GTDIM) {
        sum += correlation_function(g_in[i], g_in0[i], q_vector);
    }
    // reduced value for this thread
    s_sum[TID] = sum;
    __syncthreads();

    // compute reduced value for all threads in block
    reduce<THREADS / 2, sum_>(sum, s_sum);

    if (TID < 1) {
        // store block reduced value in global memory
        g_sum[blockIdx.x] = sum;
    }
}

template <typename vector_type>
__device__ dsfloat incoherent_scattering_function(vector_type const& r, vector_type const& r0, vector_type const& q)
{
    // accurate trigonometric function requires local memory
    return cosf((r - r0) * q);
}

template <typename vector_type,
          void (*correlation_function)(dsfloat&, dsfloat&, vector_type const&, vector_type const&),
          typename coalesced_vector_type,
          typename uncoalesced_vector_type>
__global__ void accumulate(coalesced_vector_type const* g_in, uncoalesced_vector_type const q_vector, dsfloat* g_real, dsfloat* g_imag, uint n)
{
    __shared__ dsfloat s_real[THREADS];
    __shared__ dsfloat s_imag[THREADS];

    dsfloat real = 0;
    dsfloat imag = 0;

    // load values from global device memory
    for (uint i = GTID; i < n; i += GTDIM) {
        correlation_function(real, imag, g_in[i], q_vector);
    }
    // reduced value for this thread
    s_real[TID] = real;
    s_imag[TID] = imag;
    __syncthreads();

    // compute reduced value for all threads in block
    reduce<THREADS / 2, complex_sum_>(real, imag, s_real, s_imag);

    if (TID < 1) {
        // store block reduced value in global memory
        g_real[blockIdx.x] = real;
        g_imag[blockIdx.x] = imag;
    }
}

template <typename vector_type>
__device__ void coherent_scattering_function(dsfloat& real, dsfloat& imag, vector_type const& r, vector_type const& q)
{
    float c, s;
    // accurate trigonometric function requires local memory
    sincosf(r * q, &s, &c);
    real += c;
    imag += s;
}

}}} // namespace ljgpu::cu::tcf

namespace ljgpu { namespace gpu
{

/**
 * device function wrappers
 */
cuda::function<void (float4 const*, float4 const*, uint*, dsfloat*, dsfloat*, uint)>
    tcf<3>::mean_square_displacement(cu::tcf::accumulate<cu::vector<float, 3>, cu::tcf::mean_square_displacement>);
cuda::function<void (float4 const*, float4 const*, uint*, dsfloat*, dsfloat*, uint)>
    tcf<3>::mean_quartic_displacement(cu::tcf::accumulate<cu::vector<float, 3>, cu::tcf::mean_quartic_displacement>);
cuda::function<void (float4 const*, float4 const*, uint*, dsfloat*, dsfloat*, uint)>
    tcf<3>::velocity_autocorrelation(cu::tcf::accumulate<cu::vector<float, 3>, cu::tcf::velocity_autocorrelation>);
cuda::function<void (float4 const*, float4 const*, float3 const, dsfloat*, uint)>
    tcf<3>::incoherent_scattering_function(cu::tcf::accumulate<cu::vector<float, 3>, cu::tcf::incoherent_scattering_function>);
cuda::function<void (float4 const*, float3 const, dsfloat*, dsfloat*, uint)>
    tcf<3>::coherent_scattering_function(cu::tcf::accumulate<cu::vector<float, 3>, cu::tcf::coherent_scattering_function>);

cuda::function<void (float2 const*, float2 const*, uint*, dsfloat*, dsfloat*, uint)>
    tcf<2>::mean_square_displacement(cu::tcf::accumulate<cu::vector<float, 2>, cu::tcf::mean_square_displacement>);
cuda::function<void (float2 const*, float2 const*, uint*, dsfloat*, dsfloat*, uint)>
    tcf<2>::mean_quartic_displacement(cu::tcf::accumulate<cu::vector<float, 2>, cu::tcf::mean_quartic_displacement>);
cuda::function<void (float2 const*, float2 const*, uint*, dsfloat*, dsfloat*, uint)>
    tcf<2>::velocity_autocorrelation(cu::tcf::accumulate<cu::vector<float, 2>, cu::tcf::velocity_autocorrelation>);
cuda::function<void (float2 const*, float2 const*, float2 const, dsfloat*, uint)>
    tcf<2>::incoherent_scattering_function(cu::tcf::accumulate<cu::vector<float, 2>, cu::tcf::incoherent_scattering_function>);
cuda::function<void (float2 const*, float2 const, dsfloat*, dsfloat*, uint)>
    tcf<2>::coherent_scattering_function(cu::tcf::accumulate<cu::vector<float, 2>, cu::tcf::coherent_scattering_function>);

}} // namespace ljgpu::gpu
