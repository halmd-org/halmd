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

#include <float.h>
#include "ljfluid_glue.hpp"
#include "cutil.h"
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

/** potential cutoff distance */
static __constant__ float r_cut;
/** squared cutoff distance */
static __constant__ float rr_cut;
/** cutoff energy for Lennard-Jones potential at cutoff length */
static __constant__ float en_cut;
#ifdef USE_CELL
/** number of cells per dimension */
static __constant__ unsigned int ncell;
#endif
#ifdef USE_SMOOTH_POTENTIAL
/** squared inverse potential smoothing function scale parameter */
static __constant__ float rri_smooth;
#endif


#ifdef USE_CELL
/**
 * determine cell index for a particle
 */
template <typename T>
__device__ unsigned int compute_cell(T const& r)
{
    //
    // Mapping the positional coordinates of a particle to its corresponding
    // cell index is the most delicate part of the cell lists update.
    // The safest way is to combine round-towards-zero with a successive
    // integer modulo operation, which comes with a performance penalty.
    //
    // As an efficient alternative, we transform the coordinates to the
    // half-open unit interval [0.0, 1.0) and multiply with the number
    // of cells per dimension afterwards.
    //
    const T cell = (__saturatef(r / box) * (1.f - FLT_EPSILON)) * ncell;
#ifdef DIM_3D
    return (uint(cell.z) * ncell + uint(cell.y)) * ncell + uint(cell.x);
#else
    return uint(cell.y) * ncell + uint(cell.x);
#endif
}
#endif /* USE_CELL */

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
    // periodically reduced coordinates
    r += dr;
    r -= floorf(r / box) * box;
    // periodically extended coordinates
    R += dr;
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
__device__ void compute_force(T const& r1, T const& r2, T& f, T& ff, float& en, float& virial)
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

#ifdef USE_SMOOTH_POTENTIAL
    // compute potential smoothing function
    float p = sqrtf(rr);
    float y = p - r_cut;
    float x2 = y * y * rri_smooth;
    float x4 = x2 * x2;
    float d = 1 + x4;
    float g = x4 / d;
    float h = 4 * y * x2 / (d * d * p);

    // apply smoothing function to obtain C^1 force function
    fval = g * fval - h * pot;
    // apply smoothing function to obtain C^2 potential function
    pot = g * pot;
#endif

    // virial equation sum
    virial += 0.5f * fval * rr;
    // potential energy contribution of this particle
    en += 0.5f * pot;
    // force from other particle acting on this particle
    T fi = fval * r;

    //
    // The summation of all forces acting on a particle is the most
    // critical part what concerns longtime accuracy of the simulation.
    //
    // Naively adding all forces with a single-precision operation is fine
    // with the Lennard-Jones potential using the N-squared algorithm, as
    // the force exhibits both a repulsive and an attractive part, and the
    // particles are more or less in random order. Thus, summing over all
    // forces comprises negative and positive summands in random order.
    //
    // With the WCA potential (Weeks-Chandler-Andersen, purely repulsive
    // part of the shifted Lennard-Jones potential) using the N-squared
    // algorithm, the center of mass velocity initially set to zero
    // undergoes a peculiar drift, but finally saturates at |v| ~ 5e-6.
    // Using the cell algorithm with the WCA potential however results
    // in a continuously drifting center of mass velocity, independent
    // of the chosen simulation timestep.
    //
    // The reason for this behaviour lies in the disadvantageous summing
    // order: With a purely repulsive potential, the summed forces of a
    // single neighbour cell will more or less have the same direction.
    // Thus, when adding the force sums of all neighbour cells, we add
    // huge force sums which will mostly cancel each other out in an
    // equilibrated system, giving a small and very inaccurate total
    // force due to being limited to single-precision floating-point
    // arithmetic.
    //
    // Indeed, implemeting the summation with a native-pair arithmetic,
    //
    //   Seppo Linnainmaa, Software for Doubled-Precision Floating-Point
    //   Computations, ACM Trans. Math. Softw., ACM, 1981, 7, pp.272-283,
    //
    // remedies the velocity drift, with fluctuations < 1e-7, as would
    // be expected from a symplectic integrator such as the Verlet
    // algorithm used here.
    // 

    T z = f + fi;
    T q = f - z;
    T zz = ((q + fi) + (f - (q + z))) + ff;
    f = z + zz;
    ff = (z - f) + zz;
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

#ifdef USE_CELL
/**
 * compute neighbour cell
 */
template <typename I>
__device__ unsigned int compute_neighbour_cell(I const& offset)
{
#ifdef DIM_3D
    // cell belonging to this execution block
    I cell = make_int3(blockIdx.x % ncell, (blockIdx.x / ncell) % ncell, blockIdx.x / ncell / ncell);
    // neighbour cell of this cell
    I neighbour = make_int3((cell.x + ncell + offset.x) % ncell, (cell.y + ncell + offset.y) % ncell, (cell.z + ncell + offset.z) % ncell);

    return (neighbour.z * ncell + neighbour.y) * ncell + neighbour.x;
#else
    // cell belonging to this execution block
    I cell = make_int2(blockIdx.x % ncell, blockIdx.x / ncell);
    // neighbour cell of this cell
    I neighbour = make_int2((cell.x + ncell + offset.x) % ncell, (cell.y + ncell + offset.y) % ncell);

    return neighbour.y * ncell + neighbour.x;
#endif
}

/**
 * compute forces with particles in a neighbour cell
 */
template <unsigned int block_size, bool same_cell, typename T, typename U, typename I>
__device__ void compute_cell_forces(U const* g_r, int const* g_n, I const& offset, T const& r, int const& n, T& f, T& ff, float& en, float& virial)
{
    __shared__ T s_r[block_size];
    __shared__ int s_n[block_size];

    // compute cell index
    unsigned int cell = compute_neighbour_cell(offset);

    // load particles coordinates for cell
    s_r[threadIdx.x] = unpack(g_r[cell * block_size + threadIdx.x]);
    s_n[threadIdx.x] = g_n[cell * block_size + threadIdx.x];
    __syncthreads();

    if (IS_REAL_PARTICLE(n)) {
	for (unsigned int i = 0; i < block_size; ++i) {
	    // skip placeholder particles
	    if (!IS_REAL_PARTICLE(s_n[i]))
		break;
	    // skip same particle
	    if (same_cell && threadIdx.x == i)
		continue;

	    compute_force(r, s_r[i], f, ff, en, virial);
	}
    }
    __syncthreads();
}

/**
 * n-dimensional MD simulation step
 */
template <unsigned int block_size, typename T, typename U>
__global__ void mdstep(U const* g_r, U* g_v, U* g_f, int const* g_tag, float* g_en, float* g_virial)
{
    // load particle associated with this thread
    T r = unpack(g_r[GTID]);
    T v = unpack(g_v[GTID]);
    int tag = g_tag[GTID];

    // potential energy contribution
    float en = 0;
    // virial equation sum contribution
    float virial = 0;

#ifdef DIM_3D
    T f = make_float3(0, 0, 0);
    T ff = make_float3(0, 0, 0);
#else
    T f = make_float2(0, 0);
    T ff = make_float2(0, 0);
#endif

    //
    // The summation over all forces acting on a particle is the most
    // delicate part of the MD step. As already explained above, with
    // the cell lists algorithm, summing over partial cell force sums
    // in an inappropriate order will cause a continuous drift of the
    // center of mass velocity due to single precision floating-point
    // arithmetic.
    //
    // Besides implementing the summation in double precision arithmetic,
    // choosing the order of summation over cells such that one partial
    // neighbour cell force sum is always followed by the sum of the
    // opposite neighbour cell also remedies the velocity drift.
    //

#ifdef DIM_3D
    // sum forces over this cell
    compute_cell_forces<block_size, true>(g_r, g_tag, make_int3( 0,  0,  0), r, tag, f, ff, en, virial);
    // sum forces over 26 neighbour cells, grouped into 13 pairs of mutually opposite cells
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3(-1, -1, -1), r, tag, f, ff, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3(+1, +1, +1), r, tag, f, ff, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3(-1, -1, +1), r, tag, f, ff, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3(+1, +1, -1), r, tag, f, ff, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3(-1, +1, +1), r, tag, f, ff, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3(+1, -1, -1), r, tag, f, ff, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3(+1, -1, +1), r, tag, f, ff, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3(-1, +1, -1), r, tag, f, ff, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3(-1, -1,  0), r, tag, f, ff, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3(+1, +1,  0), r, tag, f, ff, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3(-1, +1,  0), r, tag, f, ff, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3(+1, -1,  0), r, tag, f, ff, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3(-1,  0, -1), r, tag, f, ff, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3(+1,  0, +1), r, tag, f, ff, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3(-1,  0, +1), r, tag, f, ff, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3(+1,  0, -1), r, tag, f, ff, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3( 0, -1, -1), r, tag, f, ff, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3( 0, +1, +1), r, tag, f, ff, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3( 0, -1, +1), r, tag, f, ff, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3( 0, +1, -1), r, tag, f, ff, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3(-1,  0,  0), r, tag, f, ff, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3(+1,  0,  0), r, tag, f, ff, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3( 0, -1,  0), r, tag, f, ff, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3( 0, +1,  0), r, tag, f, ff, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3( 0,  0, -1), r, tag, f, ff, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3( 0,  0, +1), r, tag, f, ff, en, virial);
#else
    // sum forces over this cell
    compute_cell_forces<block_size, true>(g_r, g_tag, make_int2( 0,  0), r, tag, f, ff, en, virial);
    // sum forces over 8 neighbour cells, grouped into 4 pairs of mutually opposite cells
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int2(-1, -1), r, tag, f, ff, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int2(+1, +1), r, tag, f, ff, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int2(-1, +1), r, tag, f, ff, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int2(+1, -1), r, tag, f, ff, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int2(-1,  0), r, tag, f, ff, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int2(+1,  0), r, tag, f, ff, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int2( 0, -1), r, tag, f, ff, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int2( 0, +1), r, tag, f, ff, en, virial);
#endif

    // second leapfrog step as part of integration of equations of motion
    leapfrog_full_step(v, f);

    // store particle associated with this thread
    g_v[GTID] = pack(v);
    g_f[GTID] = pack(f);
    g_en[GTID] = en;
    g_virial[GTID] = virial;
}

#else /* USE_CELL */

/**
 * MD simulation step
 */
template <typename T, typename U>
__global__ void mdstep(U* g_r, U* g_v, U* g_f, float* g_en, float* g_virial)
{
    extern __shared__ T s_r[];

    // load particle associated with this thread
    T r = unpack(g_r[GTID]);
    T v = unpack(g_v[GTID]);

    // potential energy contribution
    float en = 0;
    // virial equation sum contribution
    float virial = 0;

#ifdef DIM_3D
    T f = make_float3(0, 0, 0);
    T ff = make_float3(0, 0, 0);
#else
    T f = make_float2(0, 0);
    T ff = make_float2(0, 0);
#endif

    // iterate over all blocks
    for (unsigned int k = 0; k < gridDim.x; k++) {
	// load positions of particles within block
	s_r[TID] = unpack(g_r[k * blockDim.x + TID]);
	__syncthreads();

	// iterate over all particles within block
	for (unsigned int j = 0; j < blockDim.x; j++) {
	    // skip placeholder particles
	    if (k * blockDim.x + j >= npart)
		continue;
	    // skip identical particle
	    if (blockIdx.x == k && TID == j)
		continue;

	    // compute Lennard-Jones force with particle
	    compute_force(r, s_r[j], f, ff, en, virial);
	}
	__syncthreads();
    }

    // second leapfrog step of integration of equations of motion
    leapfrog_full_step(v, f);

    // store particle associated with this thread
    g_v[GTID] = pack(v);
    g_f[GTID] = pack(f);
    g_en[GTID] = en;
    g_virial[GTID] = virial;
}

#endif /* USE_CELL */

/**
 * place particles on a face centered cubic (FCC) lattice
 */
template <typename T, typename U>
__global__ void lattice(U* g_r, unsigned int n)
{
    T r;
#ifdef DIM_3D
    // compose primitive vectors from 1-dimensional index
    r.x = ((GTID >> 2) % n) + ((GTID ^ (GTID >> 1)) & 1) / 2.f;
    r.y = ((GTID >> 2) / n % n) + (GTID & 1) / 2.f;
    r.z = ((GTID >> 2) / n / n) + (GTID & 2) / 4.f;
#else
    // compose primitive vectors from 1-dimensional index
    r.x = ((GTID >> 1) % n) + (GTID & 1) / 2.f;
    r.y = ((GTID >> 1) / n) + (GTID & 1) / 2.f;
#endif

    g_r[GTID] = pack(r * (box / n));
}

/**
 * generate random n-dimensional Maxwell-Boltzmann distributed velocities
 */
template <typename U>
__global__ void boltzmann(U* g_v, float temp, ushort3* g_rng)
{
    ushort3 rng = g_rng[GTID];
    U v;

    // Box-Muller transformation generates two variates at once
    rand48::gaussian(v.x, v.y, temp, rng);
#ifdef DIM_3D
    rand48::gaussian(v.z, v.w, temp, rng);
#endif

    g_rng[GTID] = rng;
    g_v[GTID] = v;
}

#ifdef USE_CELL
/**
 * assign particles to cells
 */
template <unsigned int cell_size, typename T, typename U>
__global__ void assign_cells(U const* g_part, U* g_r, int* g_tag)
{
    __shared__ T s_block[cell_size];
    __shared__ T s_r[cell_size];
    __shared__ int s_icell[cell_size];
    __shared__ int s_tag[cell_size];
    // number of particles in cell
    unsigned int n = 0;

    // mark all particles in cell as virtual particles
    s_tag[threadIdx.x] = VIRTUAL_PARTICLE;

    __syncthreads();

    for (unsigned int i = 0; i < npart; i += cell_size) {
	// load block of particles from global device memory
	T r = unpack(g_part[i + threadIdx.x]);
	s_block[threadIdx.x] = r;
	s_icell[threadIdx.x] = compute_cell(r);
	__syncthreads();

	if (threadIdx.x == 0) {
	    for (unsigned int j = 0; j < cell_size && (i + j) < npart; j++) {
		if (s_icell[j] == blockIdx.x) {
		    // store particle in cell
		    s_r[n] = s_block[j];
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
    g_r[blockIdx.x * cell_size + threadIdx.x] = pack(s_r[threadIdx.x]);
    g_tag[blockIdx.x * cell_size + threadIdx.x] = s_tag[threadIdx.x];
}

/**
 * examine neighbour cell for particles which moved into this block's cell
 */
template <unsigned int cell_size, typename T, typename U, typename I>
__device__ void examine_cell(I const& offset, U const* g_ir, U const* g_iR, U const* g_iv, int const* g_itag, T* s_or, T* s_oR, T* s_ov, int* s_otag, unsigned int& npart)
{
    __shared__ T s_ir[cell_size];
    __shared__ T s_iR[cell_size];
    __shared__ T s_iv[cell_size];
    __shared__ int s_itag[cell_size];
    __shared__ unsigned int s_cell[cell_size];

    // compute cell index
    unsigned int cell = compute_neighbour_cell(offset);

    // load particles in cell from global device memory
    T r = unpack(g_ir[cell * cell_size + threadIdx.x]);
    s_ir[threadIdx.x] = r;
    s_iR[threadIdx.x] = unpack(g_iR[cell * cell_size + threadIdx.x]);
    s_iv[threadIdx.x] = unpack(g_iv[cell * cell_size + threadIdx.x]);
    s_itag[threadIdx.x] = g_itag[cell * cell_size + threadIdx.x];
    // compute new cell
    s_cell[threadIdx.x] = compute_cell(r);
    __syncthreads();

    if (threadIdx.x == 0) {
	for (unsigned int j = 0; j < cell_size; j++) {
	    // skip virtual particles
	    if (!IS_REAL_PARTICLE(s_itag[j]))
		break;

	    // if particle belongs to this cell
	    if (s_cell[j] == blockIdx.x && npart < cell_size) {
		// store particle in cell
		s_or[npart] = s_ir[j];
		s_oR[npart] = s_iR[j];
		s_ov[npart] = s_iv[j];
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
template <unsigned int cell_size, typename T, typename U>
__global__ void update_cells(U const* g_ir, U const* g_iR, U const* g_iv, int const* g_itag, U* g_or, U* g_oR, U* g_ov, int* g_otag)
{
    __shared__ T s_or[cell_size];
    __shared__ T s_oR[cell_size];
    __shared__ T s_ov[cell_size];
    __shared__ int s_otag[cell_size];
    // number of particles in cell
    unsigned int n = 0;

    // mark all particles in cell as virtual particles
    s_otag[threadIdx.x] = VIRTUAL_PARTICLE;
    __syncthreads();

#ifdef DIM_3D
    examine_cell<cell_size>(make_int3( 0,  0,  0), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    // visit 26 neighbour cells
    examine_cell<cell_size>(make_int3(-1,  0,  0), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3(+1,  0,  0), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3( 0, -1,  0), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3( 0, +1,  0), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3(-1, -1,  0), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3(-1, +1,  0), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3(+1, -1,  0), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3(+1, +1,  0), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3( 0,  0, -1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3(-1,  0, -1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3(+1,  0, -1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3( 0, -1, -1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3( 0, +1, -1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3(-1, -1, -1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3(-1, +1, -1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3(+1, -1, -1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3(+1, +1, -1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3( 0,  0, +1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3(-1,  0, +1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3(+1,  0, +1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3( 0, -1, +1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3( 0, +1, +1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3(-1, -1, +1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3(-1, +1, +1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3(+1, -1, +1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int3(+1, +1, +1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
#else
    examine_cell<cell_size>(make_int2( 0,  0), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    // visit 8 neighbour cells
    examine_cell<cell_size>(make_int2(-1,  0), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int2(+1,  0), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int2( 0, -1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int2( 0, +1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int2(-1, -1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int2(-1, +1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int2(+1, -1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
    examine_cell<cell_size>(make_int2(+1, +1), g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);
#endif

    // store cell in global device memory
    g_or[blockIdx.x * cell_size + threadIdx.x] = pack(s_or[threadIdx.x]);
    g_oR[blockIdx.x * cell_size + threadIdx.x] = pack(s_oR[threadIdx.x]);
    g_ov[blockIdx.x * cell_size + threadIdx.x] = pack(s_ov[threadIdx.x]);
    g_otag[blockIdx.x * cell_size + threadIdx.x] = s_otag[threadIdx.x];
}
#endif  /* USE_CELL */

} // namespace mdsim


namespace mdsim { namespace gpu { namespace ljfluid
{

#ifdef DIM_3D
function<void (float4*, float4*, float4*, float4 const*)> inteq(mdsim::inteq<float3>);
#ifdef USE_CELL
function<void (float4 const*, float4*, float4*, int const*, float*, float*)> mdstep(mdsim::mdstep<CELL_SIZE, float3>);
function<void (float4 const*, float4*, int*)> assign_cells(mdsim::assign_cells<CELL_SIZE, float3>);
function<void (float4 const*, float4 const*, float4 const*, int const*, float4*, float4*, float4*, int*)> update_cells(mdsim::update_cells<CELL_SIZE, float3>);
#else
function<void (float4*, float4*, float4*, float*, float*)> mdstep(mdsim::mdstep<float3>);
#endif
function<void (float4*, unsigned int)> lattice(mdsim::lattice<float3>);
function<void (float4*, float, ushort3*)> boltzmann(mdsim::boltzmann);
#else /* DIM_3D */
function<void (float2*, float2*, float2*, float2 const*)> inteq(mdsim::inteq<float2>);
#ifdef USE_CELL
function<void (float2 const*, float2*, float2*, int const*, float*, float*)> mdstep(mdsim::mdstep<CELL_SIZE, float2>);
function<void (float2 const*, float2*, int*)> assign_cells(mdsim::assign_cells<CELL_SIZE, float2>);
function<void (float2 const*, float2 const*, float2 const*, int const*, float2*, float2*, float2*, int*)> update_cells(mdsim::update_cells<CELL_SIZE, float2>);
#else
function<void (float2*, float2*, float2*, float*, float*)> mdstep(mdsim::mdstep<float2>);
#endif
function<void (float2*, unsigned int)> lattice(mdsim::lattice<float2>);
function<void (float2*, float, ushort3*)> boltzmann(mdsim::boltzmann);
#endif /* DIM_3D */

symbol<unsigned int> npart(mdsim::npart);
symbol<float> box(mdsim::box);
symbol<float> timestep(mdsim::timestep);
symbol<float> r_cut(mdsim::r_cut);
symbol<float> rr_cut(mdsim::rr_cut);
symbol<float> en_cut(mdsim::en_cut);
#ifdef USE_CELL
symbol<unsigned int> ncell(mdsim::ncell);
#endif
#ifdef USE_SMOOTH_POTENTIAL
symbol<float> rri_smooth(mdsim::rri_smooth);
#endif

symbol<uint3> a(::rand48::a);
symbol<uint3> c(::rand48::c);

}}} // namespace mdsim::gpu::ljfluid
