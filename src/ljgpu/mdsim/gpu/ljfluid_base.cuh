/* Lennard-Jones fluid kernel
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
#include <ljgpu/math/gpu/dsfun.cuh>
#include <ljgpu/math/gpu/vector2d.cuh>
#include <ljgpu/math/gpu/vector3d.cuh>

namespace ljgpu { namespace gpu
{

/** number of particles */
static __constant__ uint npart;

/** simulation timestemp */
static __constant__ float timestep;
/** periodic box length */
static __constant__ float box;

/** potential cutoff radius */
static __constant__ float r_cut;
/** squared cutoff radius */
static __constant__ float rr_cut;
/** cutoff energy for Lennard-Jones potential at cutoff length */
static __constant__ float en_cut;

/** squared inverse potential smoothing function scale parameter */
static __constant__ float rri_smooth;

/**
 * first leapfrog step of integration of equations of motion
 */
template <typename T>
__device__ void leapfrog_half_step(T& r_, T& r, T& v, T const& f)
{
    // half step velocity
    v += f * (timestep / 2);
    // full step coordinates
    T dr = v * timestep;
    // periodically reduced coordinates
    r_ += dr;
    r_ -= floorf(r_ / box) * box;
    // periodically extended coordinates
    r += dr;
}

/**
 * second leapfrog step of integration of equations of motion
 */
template <typename T>
__device__ void leapfrog_full_step(T& v, T const& f)
{
    // full step velocity
    v += f * (timestep / 2);
}

/**
 * calculate potential smoothing function and its first derivative
 *
 * returns tuple (r, h(r), h'(r))
 */
__device__ float3 compute_smooth_function(float const& r)
{
    float y = r - r_cut;
    float x2 = y * y * rri_smooth;
    float x4 = x2 * x2;
    float x4i = 1 / (1 + x4);
    float3 h;
    h.x = r;
    h.y = x4 * x4i;
    h.z = 4 * y * x2 * x4i * x4i;
    return h;
}

/**
 * sample potential smoothing function in given range
 */
__global__ void sample_smooth_function(float3* g_h, const float2 r)
{
    g_h[GTID] = compute_smooth_function(r.x + (r.y - r.x) / GTDIM * GTID);
}

/**
 * calculate particle force using Lennard-Jones potential
 */
template <typename T, typename TT>
__device__ void compute_force(T const& r1, T const& r2, TT& f, float& en, float& virial)
{
    // particle distance vector
    T r = r1 - r2;
    // enforce periodic boundary conditions
    r -= rintf(__fdividef(r, box)) * box;
    // squared particle distance
    float rr = r * r;

    // enforce cutoff length
    if (rr >= rr_cut) return;

    // compute Lennard-Jones force in reduced units
    float rri = 1 / rr;
    float ri6 = rri * rri * rri;
    float fval = 48 * rri * ri6 * (ri6 - 0.5f);
    // compute shifted Lennard-Jones potential
    float pot = 4 * ri6 * (ri6 - 1) - en_cut;
#ifdef USE_POTENTIAL_SMOOTHING
    // compute smoothing function and its first derivative
    const float3 h = compute_smooth_function(sqrtf(rr));
    // apply smoothing function to obtain C^1 force function
    fval = h.y * fval - h.z * pot / h.x;
    // apply smoothing function to obtain C^2 potential function
    pot = h.y * pot;
#endif

    // virial equation sum
    virial += 0.5f * fval * rr;
    // potential energy contribution of this particle
    en += 0.5f * pot;
    // force from other particle acting on this particle
    f += fval * r;
}

/**
 * first leapfrog step of integration of equations of motion
 */
template <typename T, typename U>
__global__ void inteq(U* g_r, U* g_R, U* g_v, U const* g_f)
{
    T r = unpack(g_r[GTID]);
    T R = unpack(g_R[GTID]);
    T v = unpack(g_v[GTID]);
    T f = unpack(g_f[GTID]);

    leapfrog_half_step(r, R, v, f);

    g_r[GTID] = pack(r);
    g_R[GTID] = pack(R);
    g_v[GTID] = pack(v);
}

/**
 * blockwise potential energy sum
 */
__global__ void potential_energy_sum(float const* g_en, float2* g_en_sum)
{
    // single-double floating point arithmetic
    extern __shared__ dfloat s_en[];

    // load particles from global device memory
    dfloat en = 0;
    for (uint i = GTID; i < npart; i += GTDIM) {
	en += g_en[i];
    }
    // potential energy sum for this thread
    s_en[TID] = en;
    __syncthreads();

    // compute potential energy sum for all threads in block
    if (TID < 256) {
	en = en + s_en[TID + 256];
	s_en[TID] = en;
    }
    __syncthreads();
    if (TID < 128) {
	en = en + s_en[TID + 128];
	s_en[TID] = en;
    }
    __syncthreads();
    if (TID < 64) {
	en = en + s_en[TID + 64];
	s_en[TID] = en;
    }
    __syncthreads();
    if (TID < 32) {
	en = en + s_en[TID + 32];
	s_en[TID] = en;
    }
    // no further syncs needed within execution warp of 32 threads
    if (TID < 16) {
	en = en + s_en[TID + 16];
	s_en[TID] = en;
    }
    if (TID < 8) {
	en = en + s_en[TID + 8];
	s_en[TID] = en;
    }
    if (TID < 4) {
	en = en + s_en[TID + 4];
	s_en[TID] = en;
    }
    if (TID < 2) {
	en = en + s_en[TID + 2];
	s_en[TID] = en;
    }
    if (TID < 1) {
	en = en + s_en[TID + 1];
	// store potential energy block sum in global memory
	g_en_sum[blockIdx.x] = make_float2(en.f0, en.f1);
    }
}

/**
 * place particles on a face centered cubic lattice (fcc)
 */
__global__ void lattice(float4* g_r, uint n)
{
    float3 r;
    // compose primitive vectors from 1-dimensional index
    r.x = ((GTID >> 2) % n) + ((GTID ^ (GTID >> 1)) & 1) / 2.f;
    r.y = ((GTID >> 2) / n % n) + (GTID & 1) / 2.f;
    r.z = ((GTID >> 2) / n / n) + (GTID & 2) / 4.f;
    g_r[GTID] = pack(r * (box / n));
}

__global__ void lattice(float2* g_r, uint n)
{
    float2 r;
    r.x = ((GTID >> 1) % n) + (GTID & 1) / 2.f;
    r.y = ((GTID >> 1) / n) + (GTID & 1) / 2.f;
    g_r[GTID] = pack(r * (box / n));
}

/**
 * place particles on a simple cubic lattice (scc)
 */
__global__ void lattice_simple(float4* g_r, uint n)
{
    float3 r;
    r.x = (GTID % n) + 0.5f;
    r.y = (GTID / n % n) + 0.5f;
    r.z = (GTID / n / n) + 0.5f;
    g_r[GTID] = pack(r * (box / n));
}

__global__ void lattice_simple(float2* g_r, uint n)
{
    float2 r;
    r.x = (GTID % n) + 0.5f;
    r.y = (GTID / n) + 0.5f;
    g_r[GTID] = pack(r * (box / n));
}

}} // namespace ljgpu::gpu
