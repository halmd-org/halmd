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
#define CU_NAMESPACE ljfluid
#include <ljgpu/rng/gpu/rand48.cuh>

//
// compile this code only *once* per program or dynamic library,
// to avoid name conflicts of exported CUDA device symbols!
//

namespace ljgpu { namespace cu { namespace ljfluid
{

/** number of particles */
__constant__ uint npart;

/** simulation timestemp */
__constant__ float timestep;
/** periodic box length */
__constant__ float box;

/** potential cutoff radius */
__constant__ float r_cut;
/** squared cutoff radius */
__constant__ float rr_cut;
/** cutoff energy for Lennard-Jones potential at cutoff length */
__constant__ float en_cut;

/** squared inverse potential smoothing function scale parameter */
__constant__ float rri_smooth;

/**
 * first leapfrog step of integration of equations of motion
 */
template <typename T>
__device__ void leapfrog_half_step(T& r, T& R, T& v, T const& f)
{
    // half step velocity
    v += f * (timestep / 2);
    // full step coordinates
    T dr = v * timestep;
    r += dr;
    // apply periodic boundary conditions
    T dR = floorf(r / box);
    r -= dR * box;
    // periodic box traversal vector
    R += dR;
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
 * generate 4-dimensional Maxwell-Boltzmann distributed vectors
 */
__global__ void boltzmann(float4* g_v, float temperature)
{
    float4 v;
    rand48::gaussian(v, temperature);
    g_v[GTID] = v;
}

/**
 * generate 2-dimensional Maxwell-Boltzmann distributed vectors
 */
__global__ void boltzmann(float2* g_v, float temperature)
{
    float2 v;
    rand48::gaussian(v, temperature);
    g_v[GTID] = v;
}

}}} // namespace ljgpu::cu::ljfluid
