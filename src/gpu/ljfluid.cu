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

/** number of cells per dimension */
static __constant__ unsigned int ncell;


/**
 * determine cell index for a particle
 */
template <typename T>
__device__ unsigned int compute_cell(T const& r)
{
#ifdef DIM_3D
    uint3 cell = __float2uint_rz(r / (box / ncell));
    return ((cell.z % ncell) * ncell + (cell.y % ncell)) * ncell + (cell.x % ncell);
#else
    uint2 cell = __float2uint_rz(r / (box / ncell));
    return (cell.y % ncell) * ncell + (cell.x % ncell);
#endif
}


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
__device__ void leapfrog_half_step(T& r, T& rp, T& v, T const& f)
{
    // half step velocity
    v += f * (timestep / 2);
    // full step coordinates
    T dr = v * timestep;
    r += dr;
    rp += dr;
    // enforce periodic boundary conditions
    rp -= floorf(rp / box) * box;
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


/**
 * integrate equations of motion
 */
template <typename T>
__global__ void inteq(T* r, T* rp, T* v, T* f)
{
    // first leapfrog step as part of integration of equations of motion
    leapfrog_half_step(r[GTID], rp[GTID], v[GTID], f[GTID]);
}


/**
 * compute neighbour cell
 */
template <typename T>
__device__ unsigned int compute_neighbour_cell(T const& offset)
{
#ifdef DIM_3D
    // cell belonging to this execution block
    T cell = make_int3(blockIdx.x % ncell, (blockIdx.x / ncell) % ncell, blockIdx.x / ncell / ncell);
    // neighbour cell of this cell
    T neighbour = make_int3((cell.x + ncell + offset.x) % ncell, (cell.y + ncell + offset.y) % ncell, (cell.z + ncell + offset.z) % ncell);

    return (neighbour.z * ncell + neighbour.y) * ncell + neighbour.x;
#else
    // cell belonging to this execution block
    T cell = make_int2(blockIdx.x % ncell, blockIdx.x / ncell);
    // neighbour cell of this cell
    T neighbour = make_int2((cell.x + ncell + offset.x) % ncell, (cell.y + ncell + offset.y) % ncell);

    return neighbour.y * ncell + neighbour.x;
#endif
}


/**
 * compute forces with particles in a neighbour cell
 */
template <unsigned int block_size, typename T, typename U>
__device__ void compute_cell_forces(T const* g_cell, int const* g_tag, U const& offset, T const& r, T& f, int const& tag, float& en, float& virial)
{
    __shared__ T s_cell[block_size];
    __shared__ int s_tag[block_size];

    // compute cell index
    unsigned int icell = compute_neighbour_cell(offset);

    // load particles coordinates for cell
    s_cell[threadIdx.x] = g_cell[icell * block_size + threadIdx.x];
    s_tag[threadIdx.x] = g_tag[icell * block_size + threadIdx.x];
    __syncthreads();

    // check if a real particle
    if (IS_REAL_PARTICLE(tag)) {
	for (int i = 0; i < block_size; ++i) {
	    // skip same particle
	    if (blockIdx.x == icell && threadIdx.x == i) continue;

	    // skip virtual particles
	    if (!IS_REAL_PARTICLE(s_tag[i])) continue;

	    compute_force(r, s_cell[i], f, en, virial);
	}
    }
    __syncthreads();
}


/**
 * n-dimensional MD simulation step
 */
template <unsigned int block_size, typename T>
__global__ void mdstep(T const* r_, T* v_, T* f_, int const* g_tag, float2* en_)
{
    // load particle associated with this thread
    T r = r_[GTID];
    T v = v_[GTID];
    int tag = g_tag[GTID];

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

    // calculate forces for this and neighbouring cells
#ifdef DIM_3D
    compute_cell_forces<block_size>(r_, g_tag, make_int3( 0,  0,  0), r, f, tag, en, virial);
    // visit 26 neighbour cells
    compute_cell_forces<block_size>(r_, g_tag, make_int3(-1,  0,  0), r, f, tag, en, virial);
    compute_cell_forces<block_size>(r_, g_tag, make_int3(+1,  0,  0), r, f, tag, en, virial);
    compute_cell_forces<block_size>(r_, g_tag, make_int3( 0, -1,  0), r, f, tag, en, virial);
    compute_cell_forces<block_size>(r_, g_tag, make_int3( 0, +1,  0), r, f, tag, en, virial);
    compute_cell_forces<block_size>(r_, g_tag, make_int3(-1, -1,  0), r, f, tag, en, virial);
    compute_cell_forces<block_size>(r_, g_tag, make_int3(-1, +1,  0), r, f, tag, en, virial);
    compute_cell_forces<block_size>(r_, g_tag, make_int3(+1, -1,  0), r, f, tag, en, virial);
    compute_cell_forces<block_size>(r_, g_tag, make_int3(+1, +1,  0), r, f, tag, en, virial);
    compute_cell_forces<block_size>(r_, g_tag, make_int3( 0,  0, -1), r, f, tag, en, virial);
    compute_cell_forces<block_size>(r_, g_tag, make_int3(-1,  0, -1), r, f, tag, en, virial);
    compute_cell_forces<block_size>(r_, g_tag, make_int3(+1,  0, -1), r, f, tag, en, virial);
    compute_cell_forces<block_size>(r_, g_tag, make_int3( 0, -1, -1), r, f, tag, en, virial);
    compute_cell_forces<block_size>(r_, g_tag, make_int3( 0, +1, -1), r, f, tag, en, virial);
    compute_cell_forces<block_size>(r_, g_tag, make_int3(-1, -1, -1), r, f, tag, en, virial);
    compute_cell_forces<block_size>(r_, g_tag, make_int3(-1, +1, -1), r, f, tag, en, virial);
    compute_cell_forces<block_size>(r_, g_tag, make_int3(+1, -1, -1), r, f, tag, en, virial);
    compute_cell_forces<block_size>(r_, g_tag, make_int3(+1, +1, -1), r, f, tag, en, virial);
    compute_cell_forces<block_size>(r_, g_tag, make_int3( 0,  0, +1), r, f, tag, en, virial);
    compute_cell_forces<block_size>(r_, g_tag, make_int3(-1,  0, +1), r, f, tag, en, virial);
    compute_cell_forces<block_size>(r_, g_tag, make_int3(+1,  0, +1), r, f, tag, en, virial);
    compute_cell_forces<block_size>(r_, g_tag, make_int3( 0, -1, +1), r, f, tag, en, virial);
    compute_cell_forces<block_size>(r_, g_tag, make_int3( 0, +1, +1), r, f, tag, en, virial);
    compute_cell_forces<block_size>(r_, g_tag, make_int3(-1, -1, +1), r, f, tag, en, virial);
    compute_cell_forces<block_size>(r_, g_tag, make_int3(-1, +1, +1), r, f, tag, en, virial);
    compute_cell_forces<block_size>(r_, g_tag, make_int3(+1, -1, +1), r, f, tag, en, virial);
    compute_cell_forces<block_size>(r_, g_tag, make_int3(+1, +1, +1), r, f, tag, en, virial);
#else
    compute_cell_forces<block_size>(r_, g_tag, make_int2( 0,  0), r, f, tag, en, virial);
    // visit 8 neighbour cells
    compute_cell_forces<block_size>(r_, g_tag, make_int2(-1,  0), r, f, tag, en, virial);
    compute_cell_forces<block_size>(r_, g_tag, make_int2(+1,  0), r, f, tag, en, virial);
    compute_cell_forces<block_size>(r_, g_tag, make_int2( 0, -1), r, f, tag, en, virial);
    compute_cell_forces<block_size>(r_, g_tag, make_int2( 0, +1), r, f, tag, en, virial);
    compute_cell_forces<block_size>(r_, g_tag, make_int2(-1, -1), r, f, tag, en, virial);
    compute_cell_forces<block_size>(r_, g_tag, make_int2(-1, +1), r, f, tag, en, virial);
    compute_cell_forces<block_size>(r_, g_tag, make_int2(+1, -1), r, f, tag, en, virial);
    compute_cell_forces<block_size>(r_, g_tag, make_int2(+1, +1), r, f, tag, en, virial);
#endif

    // second leapfrog step as part of integration of equations of motion
    leapfrog_full_step(v, f);

    // store particle associated with this thread
    v_[GTID] = v;
    f_[GTID] = f;
    en_[GTID] = make_float2(en, virial);
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
__global__ void boltzmann(T* v_, float temp, ushort3* rng_)
{
    T v;
    ushort3 rng = rng_[GTID];

    rand48::gaussian(v.x, v.y, temp, rng);
#ifdef DIM_3D
    // Box-Muller transformation strictly generates 2 variates at once
    rand48::gaussian(v.y, v.z, temp, rng);
#endif

    rng_[GTID] = rng;
    v_[GTID] = v;
}


/**
 * assign particles to cells
 */
template <unsigned int cell_size, typename T>
__global__ void assign_cells(T const* g_part, T* g_cell, int* g_tag)
{
    __shared__ T s_block[cell_size];
    __shared__ int s_icell[cell_size];

    __shared__ T s_cell[cell_size];
    __shared__ int s_tag[cell_size];
    // number of particles in cell
    unsigned int n = 0;

    // mark all particles in cell as virtual particles
    s_tag[threadIdx.x] = VIRTUAL_PARTICLE;
    __syncthreads();

    for (unsigned int i = 0; i < npart; i += cell_size) {
	// load block of particles from global device memory
	T r = g_part[i + threadIdx.x];
	s_block[threadIdx.x] = r;
	s_icell[threadIdx.x] = compute_cell(r);
	__syncthreads();

	if (threadIdx.x == 0) {
	    for (unsigned int j = 0; j < cell_size; j++) {
		if (s_icell[j] == blockIdx.x) {
		    // store particle in cell
		    s_cell[n] = s_block[j];
		    // store particle number
		    s_tag[n] = REAL_PARTICLE(i + j);
		    // increment particle count in cell
		    ++n;
		}
	    }
	}
	__syncthreads();
    }

    // store cell in global device memory
    g_cell[blockIdx.x * cell_size + threadIdx.x] = s_cell[threadIdx.x];
    g_tag[blockIdx.x * cell_size + threadIdx.x] = s_tag[threadIdx.x];
}


/**
 * examine neighbour cell for particles which moved into this block's cell
 */
template <unsigned int cell_size, typename T, typename U>
__device__ void examine_cell(U const& offset, T const* g_ir, T const* g_irp, T const* g_iv, T const* g_if, int const* g_itag, T* s_or, T* s_orp, T* s_ov, T* s_of, int* s_otag, unsigned int& npart)
{
    __shared__ T s_ir[cell_size];
    __shared__ T s_irp[cell_size];
    __shared__ T s_iv[cell_size];
    __shared__ T s_if[cell_size];
    __shared__ int s_itag[cell_size];
    __shared__ unsigned int s_icell[cell_size];

    // compute cell index
    unsigned int icell = compute_neighbour_cell(offset);

    // load particles in cell from global device memory
    s_ir[threadIdx.x] = g_ir[icell * cell_size + threadIdx.x];
    T rp = g_irp[icell * cell_size + threadIdx.x];
    s_irp[threadIdx.x] = rp;
    s_iv[threadIdx.x] = g_iv[icell * cell_size + threadIdx.x];
    s_if[threadIdx.x] = g_if[icell * cell_size + threadIdx.x];
    s_itag[threadIdx.x] = g_itag[icell * cell_size + threadIdx.x];
    // compute new cell
    s_icell[threadIdx.x] = compute_cell(rp);
    __syncthreads();

    if (threadIdx.x == 0) {
	for (unsigned int j = 0; j < cell_size; j++) {
	    // skip virtual particles
	    if (!IS_REAL_PARTICLE(s_itag[j])) continue;

	    // if particle belongs to this cell
	    if (s_icell[j] == blockIdx.x) {
		// store particle in cell
		s_or[npart] = s_ir[j];
		s_orp[npart] = s_irp[j];
		s_ov[npart] = s_iv[j];
		s_of[npart] = s_if[j];
		s_otag[npart] = s_itag[j];
		// increment particle count in cell
		++npart;
	    }
	}
    }
    __syncthreads();
}


/**
 * update cells
 */
template <unsigned int cell_size, typename T>
__global__ void update_cells(T const* g_ir, T const* g_irp, T const* g_iv, T const* g_if, int const* g_itag, T* g_or, T* g_orp, T* g_ov, T* g_of, int* g_otag)
{
    __shared__ T s_or[cell_size];
    __shared__ T s_orp[cell_size];
    __shared__ T s_ov[cell_size];
    __shared__ T s_of[cell_size];
    __shared__ int s_otag[cell_size];
    // number of particles in cell
    unsigned int n = 0;

    // mark all particles in cell as virtual particles
    s_otag[threadIdx.x] = VIRTUAL_PARTICLE;
    __syncthreads();

#ifdef DIM_3D
    examine_cell<cell_size>(make_int3( 0,  0,  0), g_ir, g_irp, g_iv, g_if, g_itag, s_or, s_orp, s_ov, s_of, s_otag, n);
    // visit 26 neighbour cells
    examine_cell<cell_size>(make_int3(-1,  0,  0), g_ir, g_irp, g_iv, g_if, g_itag, s_or, s_orp, s_ov, s_of, s_otag, n);
    examine_cell<cell_size>(make_int3(+1,  0,  0), g_ir, g_irp, g_iv, g_if, g_itag, s_or, s_orp, s_ov, s_of, s_otag, n);
    examine_cell<cell_size>(make_int3( 0, -1,  0), g_ir, g_irp, g_iv, g_if, g_itag, s_or, s_orp, s_ov, s_of, s_otag, n);
    examine_cell<cell_size>(make_int3( 0, +1,  0), g_ir, g_irp, g_iv, g_if, g_itag, s_or, s_orp, s_ov, s_of, s_otag, n);
    examine_cell<cell_size>(make_int3(-1, -1,  0), g_ir, g_irp, g_iv, g_if, g_itag, s_or, s_orp, s_ov, s_of, s_otag, n);
    examine_cell<cell_size>(make_int3(-1, +1,  0), g_ir, g_irp, g_iv, g_if, g_itag, s_or, s_orp, s_ov, s_of, s_otag, n);
    examine_cell<cell_size>(make_int3(+1, -1,  0), g_ir, g_irp, g_iv, g_if, g_itag, s_or, s_orp, s_ov, s_of, s_otag, n);
    examine_cell<cell_size>(make_int3(+1, +1,  0), g_ir, g_irp, g_iv, g_if, g_itag, s_or, s_orp, s_ov, s_of, s_otag, n);
    examine_cell<cell_size>(make_int3( 0,  0, -1), g_ir, g_irp, g_iv, g_if, g_itag, s_or, s_orp, s_ov, s_of, s_otag, n);
    examine_cell<cell_size>(make_int3(-1,  0, -1), g_ir, g_irp, g_iv, g_if, g_itag, s_or, s_orp, s_ov, s_of, s_otag, n);
    examine_cell<cell_size>(make_int3(+1,  0, -1), g_ir, g_irp, g_iv, g_if, g_itag, s_or, s_orp, s_ov, s_of, s_otag, n);
    examine_cell<cell_size>(make_int3( 0, -1, -1), g_ir, g_irp, g_iv, g_if, g_itag, s_or, s_orp, s_ov, s_of, s_otag, n);
    examine_cell<cell_size>(make_int3( 0, +1, -1), g_ir, g_irp, g_iv, g_if, g_itag, s_or, s_orp, s_ov, s_of, s_otag, n);
    examine_cell<cell_size>(make_int3(-1, -1, -1), g_ir, g_irp, g_iv, g_if, g_itag, s_or, s_orp, s_ov, s_of, s_otag, n);
    examine_cell<cell_size>(make_int3(-1, +1, -1), g_ir, g_irp, g_iv, g_if, g_itag, s_or, s_orp, s_ov, s_of, s_otag, n);
    examine_cell<cell_size>(make_int3(+1, -1, -1), g_ir, g_irp, g_iv, g_if, g_itag, s_or, s_orp, s_ov, s_of, s_otag, n);
    examine_cell<cell_size>(make_int3(+1, +1, -1), g_ir, g_irp, g_iv, g_if, g_itag, s_or, s_orp, s_ov, s_of, s_otag, n);
    examine_cell<cell_size>(make_int3( 0,  0, +1), g_ir, g_irp, g_iv, g_if, g_itag, s_or, s_orp, s_ov, s_of, s_otag, n);
    examine_cell<cell_size>(make_int3(-1,  0, +1), g_ir, g_irp, g_iv, g_if, g_itag, s_or, s_orp, s_ov, s_of, s_otag, n);
    examine_cell<cell_size>(make_int3(+1,  0, +1), g_ir, g_irp, g_iv, g_if, g_itag, s_or, s_orp, s_ov, s_of, s_otag, n);
    examine_cell<cell_size>(make_int3( 0, -1, +1), g_ir, g_irp, g_iv, g_if, g_itag, s_or, s_orp, s_ov, s_of, s_otag, n);
    examine_cell<cell_size>(make_int3( 0, +1, +1), g_ir, g_irp, g_iv, g_if, g_itag, s_or, s_orp, s_ov, s_of, s_otag, n);
    examine_cell<cell_size>(make_int3(-1, -1, +1), g_ir, g_irp, g_iv, g_if, g_itag, s_or, s_orp, s_ov, s_of, s_otag, n);
    examine_cell<cell_size>(make_int3(-1, +1, +1), g_ir, g_irp, g_iv, g_if, g_itag, s_or, s_orp, s_ov, s_of, s_otag, n);
    examine_cell<cell_size>(make_int3(+1, -1, +1), g_ir, g_irp, g_iv, g_if, g_itag, s_or, s_orp, s_ov, s_of, s_otag, n);
    examine_cell<cell_size>(make_int3(+1, +1, +1), g_ir, g_irp, g_iv, g_if, g_itag, s_or, s_orp, s_ov, s_of, s_otag, n);
#else
    examine_cell<cell_size>(make_int2( 0,  0), g_ir, g_irp, g_iv, g_if, g_itag, s_or, s_orp, s_ov, s_of, s_otag, n);
    // visit 8 neighbour cells
    examine_cell<cell_size>(make_int2(-1,  0), g_ir, g_irp, g_iv, g_if, g_itag, s_or, s_orp, s_ov, s_of, s_otag, n);
    examine_cell<cell_size>(make_int2(+1,  0), g_ir, g_irp, g_iv, g_if, g_itag, s_or, s_orp, s_ov, s_of, s_otag, n);
    examine_cell<cell_size>(make_int2( 0, -1), g_ir, g_irp, g_iv, g_if, g_itag, s_or, s_orp, s_ov, s_of, s_otag, n);
    examine_cell<cell_size>(make_int2( 0, +1), g_ir, g_irp, g_iv, g_if, g_itag, s_or, s_orp, s_ov, s_of, s_otag, n);
    examine_cell<cell_size>(make_int2(-1, -1), g_ir, g_irp, g_iv, g_if, g_itag, s_or, s_orp, s_ov, s_of, s_otag, n);
    examine_cell<cell_size>(make_int2(-1, +1), g_ir, g_irp, g_iv, g_if, g_itag, s_or, s_orp, s_ov, s_of, s_otag, n);
    examine_cell<cell_size>(make_int2(+1, -1), g_ir, g_irp, g_iv, g_if, g_itag, s_or, s_orp, s_ov, s_of, s_otag, n);
    examine_cell<cell_size>(make_int2(+1, +1), g_ir, g_irp, g_iv, g_if, g_itag, s_or, s_orp, s_ov, s_of, s_otag, n);
#endif

    // store cell in global device memory
    g_or[blockIdx.x * cell_size + threadIdx.x] = s_or[threadIdx.x];
    g_orp[blockIdx.x * cell_size + threadIdx.x] = s_orp[threadIdx.x];
    g_ov[blockIdx.x * cell_size + threadIdx.x] = s_ov[threadIdx.x];
    g_of[blockIdx.x * cell_size + threadIdx.x] = s_of[threadIdx.x];
    g_otag[blockIdx.x * cell_size + threadIdx.x] = s_otag[threadIdx.x];
}

} // namespace mdsim


namespace mdsim { namespace gpu { namespace ljfluid
{

#ifdef DIM_3D
function<void (float3*, float3*, float3*, float3*)> inteq(mdsim::inteq);
function<void (float3*, float, ushort3*)> boltzmann(mdsim::boltzmann);
function<void (float3 const*, float3*, float3*, int const*, float2*)> mdstep(mdsim::mdstep<CELL_SIZE>);
function<void (float3*)> lattice(mdsim::lattice);
function<void (float3 const*, float3*, int*)> assign_cells(mdsim::assign_cells<CELL_SIZE>);
function<void (float3 const*, float3 const*, float3 const*, float3 const*, int const*, float3*, float3*, float3*, float3*, int*)> update_cells(mdsim::update_cells<CELL_SIZE>);
#else
function<void (float2*, float2*, float2*, float2*)> inteq(mdsim::inteq);
function<void (float2*, float, ushort3*)> boltzmann(mdsim::boltzmann);
function<void (float2 const*, float2*, float2*, int const*, float2*)> mdstep(mdsim::mdstep<CELL_SIZE>);
function<void (float2*)> lattice(mdsim::lattice);
function<void (float2 const*, float2*, int*)> assign_cells(mdsim::assign_cells<CELL_SIZE>);
function<void (float2 const*, float2 const*, float2 const*, float2 const*, int const*, float2*, float2*, float2*, float2*, int*)> update_cells(mdsim::update_cells<CELL_SIZE>);
#endif

symbol<unsigned int> npart(mdsim::npart);
symbol<float> box(mdsim::box);
symbol<float> timestep(mdsim::timestep);
symbol<float> rr_cut(mdsim::rr_cut);
symbol<float> en_cut(mdsim::en_cut);
symbol<unsigned int> ncell(mdsim::ncell);

symbol<uint3> a(::rand48::a);
symbol<uint3> c(::rand48::c);

}}} // namespace mdsim::gpu::ljfluid
