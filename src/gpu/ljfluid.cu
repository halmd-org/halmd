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

/** simulation timestemp */
static __constant__ float timestep;
/** periodic box length */
static __constant__ float box;

/** squared cutoff length */
static __constant__ float rr_cut;
/** cutoff energy for Lennard-Jones potential at cutoff length */
static __constant__ float en_cut;


/**
 * first leapfrog step of integration of equations of motion
 */
template <typename T>
__device__ void leapfrog_half_step(T& r, T& v, T f)
{
    // half step velocity
    v += f * (timestep / 2);
    // full step coordinates
    r += v * timestep;
    // enforce periodic boundary conditions
    r -= floorf(r / box) * box;
}


/**
 * second leapfrog step of integration of equations of motion
 */
template <typename T>
__device__ void leapfrog_full_step(T& v, T f)
{
    // full step velocity
    v += f * (timestep / 2);
}


/**
 * calculate particle force using Lennard-Jones potential
 */
template <typename T>
__device__ void compute_force(T r1, T r2, T& f, float& en, float& virial)
{
    // particle distance vector
    T d = r1 - r2;
    // enforce periodic boundary conditions
    d -= roundf(d / box) * box;
    // squared particle distance
    float rr = d * d;

    // enforce cutoff length
    if (rr >= rr_cut) return;

    // compute Lennard-Jones force in reduced units
    float rri = 1 / rr;
    float ri6 = rri * rri * rri;
    float fval = 48 * rri * ri6 * (ri6 - 0.5);

    // add contribution to this particle's force only
    f += fval * d;

    // potential energy contribution from this particle
    en += 2 * ri6 * (ri6 - 1) - en_cut;

    // virial equation sum
    virial += fval * rr;
}


/**
 * n-dimensional MD simulation step
 */
template <typename T>
__global__ void mdstep(T* part, T* vel, T* forces, float* en_, float* virial_)
{
    T r, v, f;
    float en, virial;

    // particles within domain
    extern __shared__ T block[];

    // load particle associated with this thread
    r = part[GTID];
    v = vel[GTID];
    f = forces[GTID];

    // first leapfrog step as part of integration of equations of motion
    leapfrog_half_step(r, v, f);

    // potential energy contribution
    en = 0.;
    // virial equation sum contribution
    virial = 0.;

    // Lennard-Jones force calculation
#ifdef DIM_3D
    f = make_float3(0., 0., 0.);
#else
    f = make_float2(0., 0.);
#endif

    for (int k = 0; k < gridDim.x; k++) {
	// load all interacting particles coordinates within domain
	block[TID] = part[k * blockDim.x + TID];

	__syncthreads();

	for (int j = 0; j < blockDim.x; j++) {
	    // same particle
	    if (blockIdx.x == k && TID == j) continue;

	    compute_force(r, block[j], f, en, virial);
	}

	__syncthreads();
    }

    // second leapfrog step as part of integration of equations of motion
    leapfrog_full_step(v, f);

    // store particle associated with this thread
    part[GTID] = r;
    vel[GTID] = v;
    forces[GTID] = f;
    en_[GTID] = en;
    virial_[GTID] = virial;
}


/**
 * place particles on an n-dimensional hypercubic lattice
 */
template <typename T>
__global__ void init_lattice(T* part, T cell, unsigned int n)
{
    T r;
    r.x = cell.x / 2 + cell.x * (GTID % n);
#ifdef DIM_3D
    r.y = cell.y / 2 + cell.y * (GTID / n % n);
    r.z = cell.z / 2 + cell.z * (GTID / n / n);
#else
    r.y = cell.y / 2 + cell.y * (GTID / n);
#endif
    part[GTID] = r;
}


/**
 * generate random n-dimensional Maxwell-Boltzmann distributed velocities
 */
template <typename T>
__global__ void init_vel(T* vel, float temp, ushort3* rng)
{
    ushort3 state = rng[GTID];
    T v;

    rand48::gaussian(v.x, v.y, temp, state);
#ifdef DIM_3D
    // Box-Muller transformation strictly generates 2 variates at once
    rand48::gaussian(v.y, v.z, temp, state);
#endif

    rng[GTID] = state;
    vel[GTID] = v;
}


/**
 * set n-dimensional force vectors to zero
 */
template <typename T>
__global__ void init_forces(T* forces)
{
#ifdef DIM_3D
    forces[GTID] = make_float3(0., 0., 0.);
#else
    forces[GTID] = make_float2(0., 0.);
#endif
}

} // namespace mdsim


namespace mdsim { namespace gpu { namespace ljfluid
{

#ifdef DIM_3D
function<void (float3*, float3*, float3*, float*, float*)> mdstep(mdsim::mdstep);
function<void (float3*, float3, unsigned int)> init_lattice(mdsim::init_lattice);
function<void (float3*, float, ushort3*)> init_vel(mdsim::init_vel);
function<void (float3*)> init_forces(mdsim::init_forces);
#else
function<void (float2*, float2*, float2*, float*, float*)> mdstep(mdsim::mdstep);
function<void (float2*, float2, unsigned int)> init_lattice(mdsim::init_lattice);
function<void (float2*, float, ushort3*)> init_vel(mdsim::init_vel);
function<void (float2*)> init_forces(mdsim::init_forces);
#endif

symbol<float> box(mdsim::box);
symbol<float> timestep(mdsim::timestep);
symbol<float> rr_cut(mdsim::rr_cut);
symbol<float> en_cut(mdsim::en_cut);

symbol<uint3> a(::rand48::a);
symbol<uint3> c(::rand48::c);

}}} // namespace mdsim::gpu::ljfluid
