/* Time correlation functions for CUDA
 *
 * Copyright © 2008-2009  Peter Colberg
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
#include <halmd/sample/gpu/tcf.hpp>

namespace halmd { namespace cu { namespace tcf
{

enum { THREADS = gpu::tcf_base::THREADS };

template <typename vector_type,
          typename correlation_function,
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
    correlation_function correlator;
    for (uint i = GTID; i < n; i += GTDIM) {
        vector_type a = g_in[i], b = g_in0[i];
        transform<accumulate_>(count, mean, variance, correlator(a, b));
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

struct mean_square_displacement
{
    template <typename vector_type>
    __device__ dsfloat operator()(vector_type const& r, vector_type const& r0) const
    {
        vector_type const dr = r - r0;
        return dr * dr;
    }
};

struct mean_quartic_displacement
{
    template <typename vector_type>
    __device__ dsfloat operator()(vector_type const& r, vector_type const& r0) const
    {
        vector_type const dr = r - r0;
        dsfloat const rr = dr * dr;
        return rr * rr;
    }
};

struct velocity_autocorrelation
{
    template <typename vector_type>
    __device__ dsfloat operator()(vector_type const& v, vector_type const& v0) const
    {
        return v * v0;
    }
};

template <typename correlation_function>
__global__ void accumulate(float const* g_in, uint* g_n, dsfloat* g_m, dsfloat* g_v, uint n)
{
    __shared__ unsigned int s_n[THREADS];
    __shared__ dsfloat s_m[THREADS];
    __shared__ dsfloat s_v[THREADS];

    unsigned int count = 0;
    dsfloat mean = 0;
    dsfloat variance = 0;

    // load values from global device memory
    correlation_function correlator;
    for (uint i = GTID; i < n; i += GTDIM) {
        if (correlator.check(i, n)) {
            dsfloat vv = g_in[i];
            transform<accumulate_>(count, mean, variance, vv);
        }
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

static __constant__ float min_fraction;

/**
 * check if particle belongs to given faction of most mobile particles
 */
struct velocity_autocorrelation_mobile
{
    __device__ bool check(unsigned int i, unsigned int npart) const
    {
        // FIXME unit test to check accumulator counts against percentage
        return i >= float2uint((__saturatef(min_fraction) * npart), cudaRoundNearest);
    }
};

static __constant__ float max_fraction;

/**
 * check if particle belongs to given faction of most immobile particles
 */
struct velocity_autocorrelation_immobile
{
    __device__ bool check(unsigned int i, unsigned int npart) const
    {
        // FIXME unit test to check accumulator counts against percentage
        return i < float2uint((__saturatef(max_fraction) * npart), cudaRoundNearest);
    }
};

template <typename vector_type,
          typename correlation_function,
          typename coalesced_vector_type,
          typename uncoalesced_vector_type>
__global__ void accumulate(coalesced_vector_type const* g_in, coalesced_vector_type const* g_in0, uncoalesced_vector_type const q_vector, dsfloat* g_sum, uint n)
{
    __shared__ dsfloat s_sum[THREADS];

    dsfloat sum = 0;

    // load values from global device memory
    correlation_function correlator;
    for (uint i = GTID; i < n; i += GTDIM) {
        vector_type a = g_in[i], b = g_in0[i], q = q_vector;
        sum += correlator(a, b, q);
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

struct incoherent_scattering_function
{
    template <typename vector_type>
    __device__ dsfloat operator()(vector_type const& r, vector_type const& r0, vector_type const& q) const
    {
        // accurate trigonometric function requires local memory
        return cosf((r - r0) * q);
    }
};

template <typename vector_type,
          typename correlation_function,
          typename coalesced_vector_type,
          typename uncoalesced_vector_type>
__global__ void accumulate(coalesced_vector_type const* g_in, uncoalesced_vector_type const q_vector, dsfloat* g_real, dsfloat* g_imag, uint n)
{
    __shared__ dsfloat s_real[THREADS];
    __shared__ dsfloat s_imag[THREADS];

    dsfloat real = 0;
    dsfloat imag = 0;

    // load values from global device memory
    correlation_function correlator;
    for (uint i = GTID; i < n; i += GTDIM) {
        vector_type a = g_in[i], q = q_vector;
        correlator(real, imag, a, q);
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

struct coherent_scattering_function
{
    template <typename vector_type>
    __device__ void operator()(dsfloat& real, dsfloat& imag, vector_type const& r, vector_type const& q) const
    {
        float c, s;
        // accurate trigonometric function requires local memory
        sincosf(r * q, &s, &c);
        real += c;
        imag += s;
    }
};

}}} // namespace halmd::cu::tcf

namespace halmd { namespace gpu
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

cuda::symbol<float>
    tcf_base::velocity_autocorrelation_mobile::min_fraction(cu::tcf::min_fraction);
cuda::symbol<float>
    tcf_base::velocity_autocorrelation_immobile::max_fraction(cu::tcf::max_fraction);
cuda::function<void (float const*, uint*, dsfloat*, dsfloat*, uint)>
    tcf_base::velocity_autocorrelation_mobile::accumulate(cu::tcf::accumulate<cu::tcf::velocity_autocorrelation_mobile>);
cuda::function<void (float const*, uint*, dsfloat*, dsfloat*, uint)>
    tcf_base::velocity_autocorrelation_immobile::accumulate(cu::tcf::accumulate<cu::tcf::velocity_autocorrelation_immobile>);

}} // namespace halmd::gpu
