/* Lennard-Jones fluid kernel
 *
 * Copyright (C) 2008  Peter Colberg
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

#include "ljfluid_glue.hpp"
#include "cutil.h"
#include "types.h"
#include "vector2d.h"
#include "vector3d.h"
#include "rand48.h"
using namespace cuda;


namespace rand48
{

/** leapfrogging multiplier */
static __constant__ uint3 a;
/** leapfrogging addend */
static __constant__ uint3 c;

}


namespace mdsim
{

/** number of particles */
static __constant__ unsigned int npart;

/** simulation timestemp */
static __constant__ float timestep;
/** periodic box length */
static __constant__ float box;

/** squared cutoff length */
static __constant__ float rr_cut;
/** cutoff energy for Lennard-Jones potential at cutoff length */
static __constant__ float en_cut;


/**
 * Verlet algorithm for integration of equations of motion
 */
template <typename T>
__device__ void verlet_step(T& r, T& rm, T& v, T const& f)
{
    T t = r;
    // update coordinates
    r = 2. * r - rm + f * (timestep * timestep);
    // update velocity
    v = (r - rm) / (2. * timestep);
    // store previous coordinates
    rm = t;
}


/**
 * first leapfrog step of integration of equations of motion
 */
template <typename T>
__device__ void leapfrog_half_step(T& r, T& v, T const& f)
{
    // half step velocity
    v += f * (timestep / 2);
    // full step coordinates
    r += v * timestep;
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
 * calculate particle force using Lennard-Jones potential
 */
template <typename T>
__device__ void compute_force(T const& r1, T const& r2, T& f, float& en, float& virial)
{
    // particle distance vector
    T r = r1 - r2;
    // enforce periodic boundary conditions
    r -= roundf(r / box) * box;
    // squared particle distance
    float rr = r * r;

    // enforce cutoff length
    if (rr >= rr_cut) return;

    // compute Lennard-Jones force in reduced units
    float rri = 1 / rr;
    float ri6 = rri * rri * rri;
    float fval = 48 * rri * ri6 * (ri6 - 0.5);

    // add contribution to this particle's force only
    f += fval * r;

    // potential energy contribution from this particle
    en += 2 * ri6 * (ri6 - 1) - en_cut;

    // virial equation sum
    virial += 0.5 * fval * rr;
}


#ifdef USE_LEAPFROG

/**
 * integrate equations of motion
 */
template <typename T>
__global__ void inteq(T* r, T* v, T* f)
{
    // first leapfrog step as part of integration of equations of motion
    leapfrog_half_step(r[GTID], v[GTID], f[GTID]);
}

#else

/**
 * integrate equations of motion
 */
template <typename T>
__global__ void inteq(T* r, T* rm, T* v, T* f)
{
    // Verlet integration of equations of motion
    verlet_step(r[GTID], rm[GTID], v[GTID], f[GTID]);
}

#endif


/**
 * n-dimensional MD simulation step
 */
template <typename T>
__global__ void mdstep(T* r_, T* v_, T* f_, float* en_, float* virial_)
{
    // particles within domain
    extern __shared__ T block[];

    // load particle associated with this thread
    T r = r_[GTID];
#ifdef USE_LEAPFROG
    T v = v_[GTID];
#endif

    // potential energy contribution
    float en = 0.;
    // virial equation sum contribution
    float virial = 0.;

    // Lennard-Jones force calculation
#ifdef DIM_3D
    T f = make_float3(0., 0., 0.);
#else
    T f = make_float2(0., 0.);
#endif

    for (int k = 0; k < gridDim.x; k++) {
	// load all interacting particles coordinates within domain
	block[TID] = r_[k * blockDim.x + TID];

	__syncthreads();

	for (int j = 0; j < blockDim.x; j++) {
	    // same particle
	    if (blockIdx.x == k && TID == j) continue;

	    compute_force(r, block[j], f, en, virial);
	}

	__syncthreads();
    }

#ifdef USE_LEAPFROG
    // second leapfrog step as part of integration of equations of motion
    leapfrog_full_step(v, f);
#endif

    // store particle associated with this thread
#ifdef USE_LEAPFROG
    v_[GTID] = v;
#endif
    f_[GTID] = f;
    en_[GTID] = en;
    virial_[GTID] = virial;
}


/**
 * place particles on a face centered cubic (FCC) lattice
 */
template <typename T>
__global__ void lattice(T* part)
{
    T r;
#ifdef DIM_3D
    // number of particles along 1 lattice dimension
    const unsigned int n = ceilf(cbrtf(npart / 4.));

    // compose primitive vectors from 1-dimensional index
    r.x = ((GTID >> 2) % n) + ((GTID ^ (GTID >> 1)) & 1) / 2.;
    r.y = ((GTID >> 2) / n % n) + (GTID & 1) / 2.;
    r.z = ((GTID >> 2) / n / n) + (GTID & 2) / 4.;
#else
    // number of particles along 1 lattice dimension
    const unsigned int n = ceilf(sqrtf(npart / 2.));

    // compose primitive vectors from 1-dimensional index
    r.x = ((GTID >> 1) % n) + (GTID & 1) / 2.;
    r.y = ((GTID >> 1) / n) + (GTID & 1) / 2.;
#endif
    part[GTID] = r * (box / n);
}


/**
 * generate random n-dimensional Maxwell-Boltzmann distributed velocities
 */
template <typename T>
#ifdef USE_LEAPFROG
__global__ void boltzmann(T* v_, float temp, ushort3* rng_)
#else
__global__ void boltzmann(T* v_, T* r_, T* rm_, float temp, ushort3* rng_)
#endif
{
#ifndef USE_LEAPFROG
    T r = r_[GTID];
#endif
    T v;
    ushort3 rng = rng_[GTID];

    rand48::gaussian(v.x, v.y, temp, rng);
#ifdef DIM_3D
    // Box-Muller transformation strictly generates 2 variates at once
    rand48::gaussian(v.y, v.z, temp, rng);
#endif

    rng_[GTID] = rng;
    v_[GTID] = v;

#ifndef USE_LEAPFROG
    // position previous time step for Verlet algorithm
    rm_[GTID] = r - v * timestep;
#endif
}

} // namespace mdsim


namespace mdsim { namespace gpu { namespace ljfluid
{

#ifdef DIM_3D

#ifdef USE_LEAPFROG
function<void (float3*, float3*, float3*)> inteq(mdsim::inteq);
function<void (float3*, float, ushort3*)> boltzmann(mdsim::boltzmann);
#else
function<void (float3*, float3*, float3*, float3*)> inteq(mdsim::inteq);
function<void (float3*, float3*, float3*, float, ushort3*)> boltzmann(mdsim::boltzmann);
#endif
function<void (float3*, float3*, float3*, float*, float*)> mdstep(mdsim::mdstep);
function<void (float3*)> lattice(mdsim::lattice);

#else

#ifdef USE_LEAPFROG
function<void (float2*, float2*, float2*)> inteq(mdsim::inteq);
function<void (float2*, float, ushort3*)> boltzmann(mdsim::boltzmann);
#else
function<void (float2*, float2*, float2*, float2*)> inteq(mdsim::inteq);
function<void (float2*, float2*, float2*, float, ushort3*)> boltzmann(mdsim::boltzmann);
#endif
function<void (float2*, float2*, float2*, float*, float*)> mdstep(mdsim::mdstep);
function<void (float2*)> lattice(mdsim::lattice);

#endif

symbol<unsigned int> npart(mdsim::npart);
symbol<float> box(mdsim::box);
symbol<float> timestep(mdsim::timestep);
symbol<float> rr_cut(mdsim::rr_cut);
symbol<float> en_cut(mdsim::en_cut);

symbol<uint3> a(::rand48::a);
symbol<uint3> c(::rand48::c);

}}} // namespace mdsim::gpu::ljfluid
