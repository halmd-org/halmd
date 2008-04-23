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
__global__ void mdstep(mdstep_param<T> state)
{
    // particles within domain
    extern __shared__ T block[];

    // load particle associated with this thread
    T r = state.r[GTID];
#ifndef USE_LEAPFROG
    T rm = state.rm[GTID];
#endif
    T v = state.v[GTID];
    T f = state.f[GTID];

#ifdef USE_LEAPFROG
    // first leapfrog step as part of integration of equations of motion
    leapfrog_half_step(r, v, f);
#else
    // Verlet integration of equations of motion
    verlet_step(r, rm, v, f);
#endif

    // potential energy contribution
    float en = 0.;
    // virial equation sum contribution
    float virial = 0.;

    // Lennard-Jones force calculation
#ifdef DIM_3D
    f = make_float3(0., 0., 0.);
#else
    f = make_float2(0., 0.);
#endif

    for (int k = 0; k < gridDim.x; k++) {
	// load all interacting particles coordinates within domain
	block[TID] = state.r[k * blockDim.x + TID];

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
    state.r[GTID] = r;
#ifndef USE_LEAPFROG
    state.rm[GTID] = rm;
#endif
    state.v[GTID] = v;
    state.f[GTID] = f;
    state.en[GTID] = en;
    state.virial[GTID] = virial;
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
__global__ void init_vel(mdstep_param<T> state, float temp, ushort3* rng)
{
#ifndef USE_LEAPFROG
    T r = state.r[GTID];
#endif
    T v;
    ushort3 rng_state = rng[GTID];

    rand48::gaussian(v.x, v.y, temp, rng_state);
#ifdef DIM_3D
    // Box-Muller transformation strictly generates 2 variates at once
    rand48::gaussian(v.y, v.z, temp, rng_state);
#endif

    rng[GTID] = rng_state;
    state.v[GTID] = v;

#ifndef USE_LEAPFROG
    // position previous time step for Verlet algorithm
    state.rm[GTID] = r - v * timestep;
#endif
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
function<void (mdstep_param<float3>)> mdstep(mdsim::mdstep);
function<void (float3*, float3, unsigned int)> init_lattice(mdsim::init_lattice);
function<void (mdstep_param<float3>, float, ushort3*)> init_vel(mdsim::init_vel);
function<void (float3*)> init_forces(mdsim::init_forces);
#else
function<void (mdstep_param<float2>)> mdstep(mdsim::mdstep);
function<void (float2*, float2, unsigned int)> init_lattice(mdsim::init_lattice);
function<void (mdstep_param<float2>, float, ushort3*)> init_vel(mdsim::init_vel);
function<void (float2*)> init_forces(mdsim::init_forces);
#endif

symbol<float> box(mdsim::box);
symbol<float> timestep(mdsim::timestep);
symbol<float> rr_cut(mdsim::rr_cut);
symbol<float> en_cut(mdsim::en_cut);

symbol<uint3> a(::rand48::a);
symbol<uint3> c(::rand48::c);

}}} // namespace mdsim::gpu::ljfluid
