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
#include <ljgpu/math/gpu/dsfloat.cuh>
#include <ljgpu/math/gpu/dsvector.cuh>
#include <ljgpu/math/gpu/vector2d.cuh>
#include <ljgpu/math/gpu/vector3d.cuh>
#include <ljgpu/math/gpu/vector4d.cuh>
#include <ljgpu/mdsim/gpu/base.hpp>
#include <ljgpu/mdsim/gpu/virial.cuh>
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
__constant__ float r_cut[3];
/** squared cutoff radius */
__constant__ float rr_cut[3];
/** Lennard-Jones potential at cutoff length in units of epsilon */
__constant__ float en_cut;
/** potential well depths in binary mixture */
__constant__ float epsilon[3];
/** squared collision diameters in binary mixture */
__constant__ float sigma2[3];

/** squared inverse potential smoothing function scale parameter */
__constant__ float rri_smooth;

/**
 * convert particle position and tag to coalesced vector type
 */
__device__ float4 attach_particle_tag(float4 r, unsigned int tag)
{
    return make_float4(r.x, r.y, r.z, __int_as_float(tag));
}

/**
 * convert coalesced vector type to particle position and tag
 */
__device__ float4 detach_particle_tag(float4 r, unsigned int& tag)
{
    tag = __float_as_int(r.w);
    return r;
}

/**
 * first leapfrog step of integration of equations of motion
 */
template <typename T, typename U, typename V>
__device__ void leapfrog_half_step(T& r, T& dr, U& R, T& v, V const& f)
{
    enum { dimension = T::static_size };
    U dR;

    // half step velocity
    v += f * (timestep / 2);
    // full step coordinates
    dr = v * timestep;
    r += dr;
    // apply periodic boundary conditions
    dR = floorf(static_cast<vector<float, dimension> >(r) / box);
    r -= dR * box;
    // periodic box traversal vector
    R += dR;
}

/**
 * second leapfrog step of integration of equations of motion
 */
template <typename T, typename U>
__device__ void leapfrog_full_step(T& v, U const& f)
{
    // full step velocity
    v += f * (timestep / 2);
}

/**
 * calculate potential smoothing function and its first derivative
 *
 * returns tuple (r, h(r), h'(r))
 */
template <mixture_type mixture, typename T>
__device__ void compute_smooth_function(T r, T& s, T& ds, unsigned int ab)
{
    T y = r - r_cut[(mixture == BINARY) ? ab : 0];
    T x2 = y * y * rri_smooth;
    T x4 = x2 * x2;
    T x4i = 1 / (1 + x4);
    // potential smoothing function
    s = x4 * x4i;
    // first derivative times (r_smooth)^(-1) [sic!]
    ds = 4 * y * rri_smooth * x2 * x4i * x4i;
}

/**
 * calculate particle force using Lennard-Jones potential
 */
template <mixture_type mixture,
	  potential_type potential,
	  typename T,
          typename U,
	  typename E,
          typename V>
__device__ void compute_force(T const& r1, T const& r2, U& f, E& en, V& virial, unsigned int ab)
{
    enum { dimension = T::static_size };
    // potential well depth
    float const eps = (mixture == BINARY) ? epsilon[ab] : 1;
    // squared collision diameter
    float const sig2 = (mixture == BINARY) ? sigma2[ab] : 1;

    // particle distance vector
    T r = r1 - r2;
    // enforce periodic boundary conditions
    r -= rintf(__fdividef(r, box)) * box;
    // squared particle distance
    typename T::value_type rr = r * r;

    // enforce cutoff length
    if (rr >= rr_cut[(mixture == BINARY) ? ab : 0]) return;

    // compute Lennard-Jones force in reduced units
    typename T::value_type rri = sig2 / rr;
    typename T::value_type ri6 = rri * rri * rri;
    typename T::value_type fval = 48 * eps * rri * ri6 * (ri6 - 0.5f) / sig2;
    // compute shifted Lennard-Jones potential
    typename T::value_type pot = (4 * ri6 * (ri6 - 1) - en_cut) * eps;

    if (potential == C2POT) {
	typename T::value_type s, ds, r_abs = sqrt(rr);
	// compute smoothing function and its first derivative
	compute_smooth_function<mixture>(r_abs, s, ds, ab);
	// apply smoothing function to obtain C¹ force function
	fval = s * fval - ds * (pot / r_abs);
	// apply smoothing function to obtain C² potential function
	pot = s * pot;
    }

    // virial equation sum
    virial += 0.5f * fval * virial::tensor(static_cast<float>(rr), static_cast<vector<float, dimension> >(r));
    // potential energy contribution of this particle
    en += 0.5f * pot;
    // force from other particle acting on this particle
    f += fval * r;
}

/**
 * sample potential smoothing function in given range
 */
__global__ void sample_smooth_function(float3* g_h, const float2 ri)
{
    float s, ds, r = ri.x + ri.y * GTID;
    compute_smooth_function<UNARY>(r, s, ds, 0);
    g_h[GTID] = make_float3(r, s, ds);
}

/**
 * sample potential and force in given range
 */
template <potential_type potential>
__global__ void sample_potential(float3* g_h, const float2 r)
{
    vector<float, 3> r1(r.x + r.y * GTID, 0, 0);
    vector<float, 3> r2(0, 0, 0);
    vector<float, 3> f = 0;
    float en = 0;
    vector<float, 4> virial = 0;
    compute_force<UNARY, potential>(r1, r2, f, en, virial, 0);
    g_h[GTID] = make_float3(r1.x, en, f.x);
}

/**
 * assign ascending particle numbers
 */
template <typename vector_type>
__global__ void init_tags(float4* g_r, unsigned int* g_tag)
{
    vector_type const r = g_r[GTID];
    unsigned int tag = VIRTUAL_PARTICLE;
    if (GTID < npart) {
	tag = GTID;
    }
    g_r[GTID] = attach_particle_tag(r, tag);
    g_tag[GTID] = tag;
}

/**
 * rescale velocities
 */
template <typename vector_type, typename T>
__global__ void rescale_velocity(T* g_v, dsfloat coeff)
{
#ifdef USE_VERLET_DSFUN
    vector<dsfloat, vector_type::static_size> v(g_v[GTID], g_v[GTID + GTDIM]);
#else
    vector_type v = g_v[GTID];
#endif
    v *= coeff;
    g_v[GTID] = static_cast<vector_type>(v);
#ifdef USE_VERLET_DSFUN
    g_v[GTID + GTDIM] = dsfloat2lo(v);
#endif
}

}}} // namespace ljgpu::cu::ljfluid

namespace ljgpu { namespace gpu
{

typedef ljfluid_base<ljfluid_impl_gpu_base> __Base;
typedef ljfluid<ljfluid_impl_gpu_base, 3> __3D;
typedef ljfluid<ljfluid_impl_gpu_base, 2> __2D;

/**
 * device constant wrappers
 */
cuda::symbol<uint> __Base::npart(cu::ljfluid::npart);
cuda::symbol<uint[]> __Base::mpart(cu::ljfluid::mpart);
cuda::symbol<float> __Base::box(cu::ljfluid::box);
cuda::symbol<float> __Base::timestep(cu::ljfluid::timestep);
cuda::symbol<float[]> __Base::r_cut(cu::ljfluid::r_cut);
cuda::symbol<float[]> __Base::rr_cut(cu::ljfluid::rr_cut);
cuda::symbol<float> __Base::en_cut(cu::ljfluid::en_cut);
cuda::symbol<float[]> __Base::epsilon(cu::ljfluid::epsilon);
cuda::symbol<float[]> __Base::sigma2(cu::ljfluid::sigma2);
cuda::symbol<float> __Base::rri_smooth(cu::ljfluid::rri_smooth);

/**
 * device function wrappers
 */
cuda::function<void (float3*, const float2)>
    __Base::sample_smooth_function(cu::ljfluid::sample_smooth_function);
cuda::function<void (float3*, const float2)>
    __Base::sample_potential(cu::ljfluid::sample_potential<C0POT>);
cuda::function<void (float3*, const float2)>
    __Base::sample_smooth_potential(cu::ljfluid::sample_potential<C2POT>);

cuda::function<void (float4*, unsigned int*)>
    __3D::init_tags(cu::ljfluid::init_tags<cu::vector<float, 3> >);
cuda::function<void (float4*, dsfloat)>
    __3D::rescale_velocity(cu::ljfluid::rescale_velocity<cu::vector<float, 3> >);

cuda::function<void (float4*, unsigned int*)>
    __2D::init_tags(cu::ljfluid::init_tags<cu::vector<float, 2> >);
cuda::function<void (float2*, dsfloat)>
    __2D::rescale_velocity(cu::ljfluid::rescale_velocity<cu::vector<float, 2> >);

}} // namespace ljgpu::gpu
