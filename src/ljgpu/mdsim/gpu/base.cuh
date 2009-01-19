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
#include <ljgpu/math/gpu/pack.cuh>
#include <ljgpu/math/gpu/vector2d.cuh>
#include <ljgpu/math/gpu/vector3d.cuh>
#include <ljgpu/mdsim/gpu/base.hpp>
#define CU_NAMESPACE ljfluid
#include <ljgpu/rng/gpu/rand48.cuh>
using namespace ljgpu::gpu;

//
// compile this code only *once* per program or dynamic library,
// to avoid name conflicts of exported CUDA device symbols!
//

namespace ljgpu { namespace cu { namespace ljfluid
{

/** number of particles */
__constant__ uint npart;
/** number of A and B particles in a binary mixture */
__constant__ uint mpart[2];

/** simulation timestemp */
__constant__ float timestep;
/** periodic box length */
__constant__ float box;

/** potential cutoff radius */
__constant__ float r_cut;
/** squared cutoff radius */
__constant__ float rr_cut[3];
/** cutoff energy for Lennard-Jones potential at cutoff length */
__constant__ float en_cut;
/** potential well depths in binary mixture */
__constant__ float epsilon[3];
/** squared collision diameters in binary mixture */
__constant__ float sigma2[3];

/** squared inverse potential smoothing function scale parameter */
__constant__ float rri_smooth;

/** heat bath coupling constant */
__constant__ float thermostat_nu;
/** heat bath temperature */
__constant__ float thermostat_temp;

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
 * random collision with heat bath
 */
template <typename T>
__device__ void anderson_thermostat(T& v)
{
    if (rand48::uniform() < (thermostat_nu * timestep)) {
	rand48::gaussian(v, thermostat_temp);
    }
}

/**
 * calculate potential smoothing function and its first derivative
 *
 * returns tuple (r, h(r), h'(r))
 */
__device__ float3 compute_smooth_function(float const& r)
{
    float y = r - /*FIXME binary mixture */ r_cut;
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
template <mixture_type mixture,
	  potential_type potential,
	  typename T,
          typename U>
__device__ void compute_force(T const& r1, T const& r2,
			      U& f, float& en, float& virial,
			      unsigned int ab = 0)
{
    // potential well depth
    float const eps = (mixture == BINARY) ? epsilon[ab] : 1;
    // squared collision diameter
    float const sig2 = (mixture == BINARY) ? sigma2[ab] : 1;

    // particle distance vector
    T r = r1 - r2;
    // enforce periodic boundary conditions
    r -= rintf(__fdividef(r, box)) * box;
    // squared particle distance
    float rr = r * r;

    // enforce cutoff length
    if (rr >= rr_cut[(mixture == BINARY) ? ab : 0]) return;

    // compute Lennard-Jones force in reduced units
    float rri = sig2 / rr;
    float ri6 = rri * rri * rri;
    float fval = 48 * eps * rri * ri6 * (ri6 - 0.5f) / sig2;
    // compute shifted Lennard-Jones potential
    float pot = (4 * ri6 * (ri6 - 1) - en_cut) * eps;

    if (potential == C2POT) {
	// compute smoothing function and its first derivative
	const float3 h = compute_smooth_function(sqrtf(rr));
	// apply smoothing function to obtain C¹ force function
	fval = h.y * fval - h.z * pot / h.x;
	// apply smoothing function to obtain C² potential function
	pot = h.y * pot;
    }

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
template <int dimension, typename T>
__global__ void inteq(float4* g_r, T* g_R, T* g_v, T const* g_f)
{
    vector<float, dimension> r, R, v, f;
    int tag;
    (r, tag) = g_r[GTID];
    R = g_R[GTID];
    v = g_v[GTID];
    f = g_f[GTID];

    leapfrog_half_step(r, R, v, f);

    g_r[GTID] = (r, tag);
    g_R[GTID] = R;
    g_v[GTID] = v;
}

/**
 * generate n-dimensional Maxwell-Boltzmann distributed vectors
 */
template <typename T>
__global__ void boltzmann(T* g_v, float temperature)
{
    T v;
    rand48::gaussian(v, temperature);
    g_v[GTID] = v;
}

/**
 * assign ascending particle numbers
 */
template <typename vector_type>
__global__ void init_tags(float4* g_r, int* g_tag)
{
    vector_type const r = g_r[GTID];
    int tag = VIRTUAL_PARTICLE;
    if (GTID < npart) {
	tag = particle_tag(GTID);
    }
    g_r[GTID] = (r, tag);
    g_tag[GTID] = tag;
}

/**
 * assign particles types in a binary mixture
 */
template <typename vector_type>
__global__ void init_types(float4* g_r, int* g_tag)
{
    vector_type const r = g_r[GTID];
    int tag = VIRTUAL_PARTICLE;
    if (GTID < mpart[0]) {
	tag = particle_tag(GTID, 0);
    }
    else if (GTID - mpart[0] < mpart[1]) {
	tag = particle_tag(GTID - mpart[0], 1);
    }
    g_r[GTID] = (r, tag);
    g_tag[GTID] = tag;
}

}}} // namespace ljgpu::cu::ljfluid

namespace ljgpu { namespace gpu
{

typedef ljfluid_base<ljfluid_impl_gpu_base> __Base;
typedef ljfluid<ljfluid_impl_gpu_base<3> > __3D;
typedef ljfluid<ljfluid_impl_gpu_base<2> > __2D;

/**
 * device constant wrappers
 */
cuda::symbol<uint> __Base::npart(cu::ljfluid::npart);
cuda::symbol<uint[]> __Base::mpart(cu::ljfluid::mpart);
cuda::symbol<float> __Base::box(cu::ljfluid::box);
cuda::symbol<float> __Base::timestep(cu::ljfluid::timestep);
cuda::symbol<float> __Base::r_cut(cu::ljfluid::r_cut);
cuda::symbol<float[]> __Base::rr_cut(cu::ljfluid::rr_cut);
cuda::symbol<float> __Base::en_cut(cu::ljfluid::en_cut);
cuda::symbol<float[]> __Base::epsilon(cu::ljfluid::epsilon);
cuda::symbol<float[]> __Base::sigma2(cu::ljfluid::sigma2);
cuda::symbol<float> __Base::rri_smooth(cu::ljfluid::rri_smooth);
cuda::symbol<float> __Base::thermostat_nu(cu::ljfluid::thermostat_nu);
cuda::symbol<float> __Base::thermostat_temp(cu::ljfluid::thermostat_temp);

cuda::symbol<uint48> __Base::rand48::a(cu::ljfluid::rand48::a);
cuda::symbol<uint48> __Base::rand48::c(cu::ljfluid::rand48::c);
cuda::symbol<ushort3*> __Base::rand48::state(cu::ljfluid::rand48::g_state);

/**
 * device function wrappers
 */
cuda::function<void (float3*, const float2)>
    __Base::sample_smooth_function(cu::ljfluid::sample_smooth_function);

cuda::function<void (float4*, float4*, float4*, float4 const*)>
    __3D::inteq(cu::ljfluid::inteq<3>);
cuda::function<void (float4*, float)>
    __3D::boltzmann(cu::ljfluid::boltzmann);
cuda::function<void (float4*, int*)>
    __3D::init_tags(cu::ljfluid::init_tags<cu::vector<float, 3> >);
cuda::function<void (float4*, int*)>
    __3D::init_types(cu::ljfluid::init_types<cu::vector<float, 3> >);

cuda::function<void (float4*, float2*, float2*, float2 const*)>
    __2D::inteq(cu::ljfluid::inteq<2>);
cuda::function<void (float2*, float)>
    __2D::boltzmann(cu::ljfluid::boltzmann);
cuda::function<void (float4*, int*)>
    __2D::init_tags(cu::ljfluid::init_tags<cu::vector<float, 2> >);
cuda::function<void (float4*, int*)>
    __2D::init_types(cu::ljfluid::init_types<cu::vector<float, 2> >);

}} // namespace ljgpu::gpu
