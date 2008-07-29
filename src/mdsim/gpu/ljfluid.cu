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
#include "algorithm.h"
#include "cutil.h"
#include "dsfun.h"
#include "vector2d.h"
#include "vector3d.h"
#include "rand48.h"


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
/** neighbour list length */
static __constant__ unsigned int nnbl;
/** potential cutoff distance with cell skin */
static __constant__ float r_cell;
/** squared potential cutoff distance with cell skin */
static __constant__ float rr_cell;
# ifdef DIM_3D
/** texture reference to periodic particle positions */
struct texture<float4, 1, cudaReadModeElementType> t_r;
# else
struct texture<float2, 1, cudaReadModeElementType> t_r;
# endif
/** number of Hilbert space-filling curve code levels */
static __constant__ unsigned int sfc_level;
#endif

#ifdef USE_SMOOTH_POTENTIAL
/** squared inverse potential smoothing function scale parameter */
static __constant__ float rri_smooth;
#endif


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

#ifdef USE_SMOOTH_POTENTIAL

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

#endif /* USE_SMOOTH_POTENTIAL */

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
#ifdef USE_SMOOTH_POTENTIAL
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

#ifdef USE_CELL

/**
 * n-dimensional MD simulation step
 */
template <typename T, typename U>
__global__ void mdstep(U const* g_r, U* g_v, U* g_f, int const* g_nbl, float* g_en, float* g_virial)
{
    // load particle associated with this thread
    const T r = unpack(g_r[GTID]);
    T v = unpack(g_v[GTID]);

    // potential energy contribution
    float en = 0;
    // virial equation sum contribution
    float virial = 0;

#ifdef DIM_3D
    dfloat3 f(make_float3(0, 0, 0));
#else
    dfloat2 f(make_float2(0, 0));
#endif

    for (unsigned int i = 0; i < nnbl; ++i) {
	// load particle number from neighbour list
	const int n = g_nbl[GTID * nnbl + i];
	// skip placeholder particles
	if (!IS_REAL_PARTICLE(n))
	    break;
	// accumulate force between particles
	compute_force(r, unpack(tex1Dfetch(t_r, n)), f, en, virial);
    }

    // second leapfrog step as part of integration of equations of motion
    leapfrog_full_step(v, f.f0);

    // store particle associated with this thread
    g_v[GTID] = pack(v);
    g_f[GTID] = pack(f.f0);
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
    dfloat3 f(make_float3(0, 0, 0));
#else
    dfloat2 f(make_float2(0, 0));
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
	    compute_force(r, s_r[j], f, en, virial);
	}
	__syncthreads();
    }

    // second leapfrog step of integration of equations of motion
    leapfrog_full_step(v, f.f0);

    // store particle associated with this thread
    g_v[GTID] = pack(v);
    g_f[GTID] = pack(f.f0);
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
 * place particles on a simple cubic (SCC) lattice
 */
template <typename T, typename U>
__global__ void lattice_simple(U* g_r, unsigned int n)
{
    T r;
#ifdef DIM_3D
    r.x = (GTID % n) + 0.5f;
    r.y = (GTID / n % n) + 0.5f;
    r.z = (GTID / n / n) + 0.5f;
#else
    r.x = (GTID % n) + 0.5f;
    r.y = (GTID / n) + 0.5f;
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
 * compute neighbour cell
 */
template <typename I>
__device__ unsigned int compute_neighbour_cell(I const& offset)
{
#ifdef DIM_3D
    const int3 cell = make_int3(blockIdx.x % ncell, (blockIdx.x / ncell) % ncell, blockIdx.x / ncell / ncell);
    const int3 neighbour = make_int3((cell.x + ncell + offset.x) % ncell, (cell.y + ncell + offset.y) % ncell, (cell.z + ncell + offset.z) % ncell);
    return (neighbour.z * ncell + neighbour.y) * ncell + neighbour.x;
#else
    const int2 cell = make_int2(blockIdx.x % ncell, blockIdx.x / ncell);
    const int2 neighbour = make_int2((cell.x + ncell + offset.x) % ncell, (cell.y + ncell + offset.y) % ncell);
    return neighbour.y * ncell + neighbour.x;
#endif
}

/**
 * update neighbour list with particles of given cell
 */
template <unsigned int cell_size, bool same_cell, typename T, typename I>
__device__ void update_cell_neighbours(I const& offset, int const* g_cell, int* g_nbl, T const& r, int const& n, unsigned int& count)
{
    __shared__ int s_n[cell_size];
    __shared__ T s_r[cell_size];

    // compute cell index
    const unsigned int cell = compute_neighbour_cell(offset);

    // load particles in cell
    s_n[threadIdx.x] = g_cell[cell * cell_size + threadIdx.x];
    s_r[threadIdx.x] = unpack(tex1Dfetch(t_r, s_n[threadIdx.x]));
    __syncthreads();

    if (IS_REAL_PARTICLE(n)) {
	for (unsigned int i = 0; i < cell_size; ++i) {
	    // particle number of cell placeholder
	    const int m = s_n[i];
	    // skip placeholder particles
	    if (!IS_REAL_PARTICLE(m))
		break;
	    // skip same particle
	    if (same_cell && i == threadIdx.x)
		continue;

	    // particle distance vector
	    T dr = r - s_r[i];
	    // enforce periodic boundary conditions
	    dr -= rintf(__fdividef(dr, box)) * box;
	    // squared particle distance
	    const float rr = dr * dr;

	    // enforce cutoff length with neighbour list skin
	    if (rr <= rr_cell && count < nnbl) {
		// uncoalesced write to neighbour list
		g_nbl[n * nnbl + count] = m;
		// increment neighbour list particle count
		count++;
	    }
	}
    }
}

/**
 * update neighbour lists
 */
template <unsigned int cell_size, typename T>
__global__ void update_neighbours(int const* g_cell, int* g_nbl)
{
    // load particle from cell placeholder
    const int n = g_cell[GTID];
    const T r = unpack(tex1Dfetch(t_r, n));
    // number of particles in neighbour list
    unsigned int count = 0;

    //
    // The summation of all forces acting on a particle is the most
    // critical part of the simulation concerning longtime accuracy.
    //
    // Naively adding all forces with a single-precision operation is fine
    // with the Lennard-Jones potential using the N-squared algorithm, as
    // the force exhibits both a repulsive and an attractive part, and the
    // particles are more or less in random order. Thus, summing over all
    // forces comprises negative and positive summands in random order.
    //
    // With the WCA potential (Weeks-Chandler-Andersen, purely repulsive
    // part of the shifted Lennard-Jones potential) using the N-squared
    // algorithm, the center of mass velocity effectively stays zero if
    // the initial list of particles arranged on a lattice is randomly
    // permuted before simulation.
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
    // Besides implementing the summation in double precision arithmetic,
    // choosing the order of summation over cells such that one partial
    // neighbour cell force sum is always followed by the sum of the
    // opposite neighbour cell softens the velocity drift.
    //

#ifdef DIM_3D
    // visit this cell
    update_cell_neighbours<cell_size, true>(make_int3( 0,  0,  0), g_cell, g_nbl, r, n, count);
    // visit 26 neighbour cells, grouped into 13 pairs of mutually opposite cells
    update_cell_neighbours<cell_size, false>(make_int3(-1, -1, -1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, false>(make_int3(+1, +1, +1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, false>(make_int3(-1, -1, +1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, false>(make_int3(+1, +1, -1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, false>(make_int3(-1, +1, +1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, false>(make_int3(+1, -1, -1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, false>(make_int3(+1, -1, +1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, false>(make_int3(-1, +1, -1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, false>(make_int3(-1, -1,  0), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, false>(make_int3(+1, +1,  0), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, false>(make_int3(-1, +1,  0), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, false>(make_int3(+1, -1,  0), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, false>(make_int3(-1,  0, -1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, false>(make_int3(+1,  0, +1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, false>(make_int3(-1,  0, +1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, false>(make_int3(+1,  0, -1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, false>(make_int3( 0, -1, -1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, false>(make_int3( 0, +1, +1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, false>(make_int3( 0, -1, +1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, false>(make_int3( 0, +1, -1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, false>(make_int3(-1,  0,  0), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, false>(make_int3(+1,  0,  0), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, false>(make_int3( 0, -1,  0), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, false>(make_int3( 0, +1,  0), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, false>(make_int3( 0,  0, -1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, false>(make_int3( 0,  0, +1), g_cell, g_nbl, r, n, count);
#else
    // visit this cell
    update_cell_neighbours<cell_size, true>(make_int2( 0,  0), g_cell, g_nbl, r, n, count);
    // visit 8 neighbour cells, grouped into 4 pairs of mutually opposite cells
    update_cell_neighbours<cell_size, false>(make_int2(-1, -1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, false>(make_int2(+1, +1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, false>(make_int2(-1, +1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, false>(make_int2(+1, -1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, false>(make_int2(-1,  0), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, false>(make_int2(+1,  0), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, false>(make_int2( 0, -1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, false>(make_int2( 0, +1), g_cell, g_nbl, r, n, count);
#endif

    // terminate neighbour list
    if (IS_REAL_PARTICLE(n) && count < nnbl) {
	g_nbl[n * nnbl + count] = VIRTUAL_PARTICLE;
    }
}

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
    return uint(cell.x) + ncell * (uint(cell.y) + ncell * uint(cell.z));
#else
    return uint(cell.x) + ncell * uint(cell.y);
#endif
}

/**
 * assign particles to cells
 */
template <unsigned int cell_size, typename U>
__global__ void assign_cells(U const* g_part, int* g_tag)
{
    __shared__ unsigned int s_cell[cell_size];
    __shared__ int s_tag[cell_size];
    // number of particles in cell
    unsigned int n = 0;

    // mark all particles in cell as virtual particles
    s_tag[threadIdx.x] = VIRTUAL_PARTICLE;
    __syncthreads();

    for (unsigned int i = 0; i < npart; i += cell_size) {
	// load block of particles from global device memory
	s_cell[threadIdx.x] = compute_cell(unpack(g_part[i + threadIdx.x]));
	__syncthreads();

	if (threadIdx.x == 0) {
	    for (unsigned int j = 0; j < cell_size && (i + j) < npart; j++) {
		if (s_cell[j] == blockIdx.x) {
		    // store particle in cell
		    s_tag[n] = REAL_PARTICLE(i + j);
		    // increment particle count in cell
		    ++n;
		}
	    }
	}
	__syncthreads();
    }

    // store cell in global device memory
    g_tag[blockIdx.x * cell_size + threadIdx.x] = s_tag[threadIdx.x];
}

/**
 * examine neighbour cell for particles which moved into this block's cell
 */
template <unsigned int cell_size, typename T>
__device__ void examine_cell(T const& offset, int const* g_itag, int* s_otag, unsigned int& npart)
{
    __shared__ int s_itag[cell_size];
    __shared__ unsigned int s_cell[cell_size];

    // compute cell index
    const unsigned int cell = compute_neighbour_cell(offset);

    // load particle numbers from global device memory
    int n = g_itag[cell * cell_size + threadIdx.x];
    s_itag[threadIdx.x] = n;

    if (IS_REAL_PARTICLE(n)) {
	// compute new cell
	s_cell[threadIdx.x] = compute_cell(unpack(tex1Dfetch(t_r, n)));
    }
    __syncthreads();

    if (threadIdx.x == 0) {
	for (unsigned int j = 0; j < cell_size; j++) {
	    // skip virtual particles
	    if (!IS_REAL_PARTICLE(s_itag[j]))
		break;

	    // if particle belongs to this cell
	    if (s_cell[j] == blockIdx.x && npart < cell_size) {
		// store particle in cell
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
template <unsigned int cell_size>
__global__ void update_cells(int const* g_itag, int* g_otag)
{
    __shared__ int s_otag[cell_size];
    // number of particles in cell
    unsigned int n = 0;

    // mark all particles in cell as virtual particles
    s_otag[threadIdx.x] = VIRTUAL_PARTICLE;
    __syncthreads();

#ifdef DIM_3D
    // visit this cell
    examine_cell<cell_size>(make_int3( 0,  0,  0), g_itag, s_otag, n);
    // visit 26 neighbour cells
    examine_cell<cell_size>(make_int3(-1, -1, -1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3(+1, +1, +1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3(-1, -1, +1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3(+1, +1, -1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3(-1, +1, +1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3(+1, -1, -1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3(+1, -1, +1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3(-1, +1, -1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3(-1, -1,  0), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3(+1, +1,  0), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3(-1, +1,  0), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3(+1, -1,  0), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3(-1,  0, -1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3(+1,  0, +1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3(-1,  0, +1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3(+1,  0, -1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3( 0, -1, -1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3( 0, +1, +1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3( 0, -1, +1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3( 0, +1, -1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3(-1,  0,  0), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3(+1,  0,  0), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3( 0, -1,  0), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3( 0, +1,  0), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3( 0,  0, -1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3( 0,  0, +1), g_itag, s_otag, n);
#else
    // visit this cell
    examine_cell<cell_size>(make_int2( 0,  0), g_itag, s_otag, n);
    // visit 8 neighbour cells
    examine_cell<cell_size>(make_int2(-1, -1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int2(+1, +1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int2(-1, +1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int2(+1, -1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int2(-1,  0), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int2(+1,  0), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int2( 0, -1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int2( 0, +1), g_itag, s_otag, n);
#endif

    // store cell in global device memory
    g_otag[blockIdx.x * cell_size + threadIdx.x] = s_otag[threadIdx.x];
}

/**
 * map n-dimensional point to 1-dimensional point on Hilbert space curve
 */
template <typename T, typename U>
__global__ void sfc_hilbert_encode(U const* g_r, unsigned int* g_sfc)
{
    //
    // Jun Wang & Jie Shan, Space-Filling Curve Based Point Clouds Index,
    // GeoComputation, 2005
    //

    // Hilbert code for particle
    unsigned int hcode = 0;

    //
    // We need to avoid ambiguities during the assignment of a particle
    // to a subcell, i.e. the particle position should never lie on an
    // edge or corner of multiple subcells, or the algorithm will have
    // trouble converging to a definite Hilbert curve.
    //
    // Therefore, we use a simple cubic lattice of predefined dimensions
    // according to the number of cells at the deepest recursion level,
    // and round the particle position to the nearest center of a cell.
    //

    // Hilbert cells per dimension at deepest recursion level
    const uint n = 1UL << sfc_level;
    // fractional index of particle's Hilbert cell in [0, n)
    const T cell = (__saturatef(unpack(g_r[GTID]) / box) * (1.f - FLT_EPSILON)) * n;

#ifdef DIM_3D

    // round particle position to center of cell in unit coordinates
    T r = make_float3(uint(cell.x) + 0.5f, uint(cell.y) + 0.5f, uint(cell.z) + 0.5f) / n;
    // use symmetric coordinates
    r -= make_float3(0.5f, 0.5f, 0.5f);

    // Hilbert vertex-to-code lookup table
    uint vv = 0x5AFA5;			// 000 001 011 010 111 110 100 101
    // Hilbert code-to-vertex lookup table
    uint a = 21;
    uint b = 18;
    uint c = 12;
    uint d = 15;
    uint e = 3;
    uint f = 0;
    uint g = 6;
    uint h = 9;

#define MASK ((1 << 3) - 1)

    // 32-bit integer for 3D Hilbert code allows a maximum of 10 levels
    for (unsigned int i = 0; i < sfc_level; ++i) {
	// determine Hilbert vertex closest to particle
	const uint x = __signbitf(r.x) & 1;
	const uint y = __signbitf(r.y) & 1;
	const uint z = __signbitf(r.z) & 1;
	// lookup Hilbert code 
	const uint v = (vv >> (3 * (x + (y << 1) + (z << 2))) & MASK);

	// scale particle coordinates to subcell
	r = 2 * r - make_float3(0.5f - x, 0.5f - y, 0.5f - z);
	// apply permutation rule according to Hilbert code
	if (v == 0) {
	    vertex_swap(vv, b, h, MASK);
	    vertex_swap(vv, c, e, MASK);
	}
	else if (v == 1 || v == 2) {
	    vertex_swap(vv, c, g, MASK);
	    vertex_swap(vv, d, h, MASK);
	}
	else if (v == 3 || v == 4) {
	    vertex_swap(vv, a, c, MASK);
#ifdef USE_ALTERNATIVE_HILBERT_3D
	    vertex_swap(vv, b, d, MASK);
	    vertex_swap(vv, e, g, MASK);
#endif
	    vertex_swap(vv, f, h, MASK);
	}
	else if (v == 5 || v == 6) {
	    vertex_swap(vv, a, e, MASK);
	    vertex_swap(vv, b, f, MASK);
	}
	else if (v == 7) {
	    vertex_swap(vv, a, g, MASK);
	    vertex_swap(vv, d, f, MASK);
	}

	// add vertex code to partial Hilbert code
	hcode = (hcode << 3) + v;
    }

#else /* ! DIM_3D */

    // round particle position to center of cell in unit coordinates
    T r = make_float2(uint(cell.x) + 0.5f, uint(cell.y) + 0.5f) / n;
    // use symmetric coordinates
    r -= make_float2(0.5f, 0.5f);

    // Hilbert vertex-to-code lookup table
    uint vv = 0x1E;			// 00 01 11 10
    // Hilbert code-to-vertex lookup table
    uint a = 6;
    uint b = 4;
    uint c = 0;
    uint d = 2;

#define MASK ((1 << 2) - 1)

    // 32-bit integer for 2D Hilbert code allows a maximum of 16 levels
    for (unsigned int i = 0; i < sfc_level; ++i) {
	// determine Hilbert vertex closest to particle
	const uint x = __signbitf(r.x) & 1;
	const uint y = __signbitf(r.y) & 1;
	// lookup Hilbert code 
	const uint v = (vv >> (2 * (x + (y << 1))) & MASK);

	// scale particle coordinates to subcell
	r = 2 * r - make_float2(0.5f - x, 0.5f - y);
	// apply permutation rule according to Hilbert code
	if (v == 0) {
	    vertex_swap(vv, b, d, MASK);
	}
	else if (v == 3) {
	    vertex_swap(vv, a, c, MASK);
	}

	// add vertex code to partial Hilbert code
	hcode = (hcode << 2) + v;
    }

#endif /* DIM_3D */
#undef MASK

    // store Hilbert code for particle
    g_sfc[GTID] = hcode;
}

/**
 * swap Hilbert spacing-filling curve vertices
 */
__device__ void vertex_swap(uint& v, uint& a, uint& b, uint const& mask)
{
    // swap bits comprising Hilbert codes in vertex-to-code lookup table
    v = (v & ~(mask << a) & ~(mask << b)) + (((v >> a) & mask) << b) + (((v >> b) & mask) << a);
    // update code-to-vertex lookup table
    swap(a, b);
}

#endif  /* USE_CELL */

} // namespace mdsim


namespace mdsim { namespace gpu { namespace ljfluid
{

#ifdef DIM_3D
cuda::function<void (float4*, float4*, float4*, float4 const*)> inteq(mdsim::inteq<float3>);
# ifdef USE_CELL
cuda::function<void (float4 const*, float4*, float4*, int const*, float*, float*)> mdstep(mdsim::mdstep<float3>);
cuda::function<void (float4 const*, int*)> assign_cells(mdsim::assign_cells<CELL_SIZE>);
cuda::function<void (int const*, int*)> update_neighbours(mdsim::update_neighbours<CELL_SIZE, float3>);
cuda::texture<float4> r(mdsim::t_r);
cuda::function<void (float4 const*, unsigned int*)> sfc_hilbert_encode(mdsim::sfc_hilbert_encode<float3>);
# else
cuda::function<void (float4*, float4*, float4*, float*, float*)> mdstep(mdsim::mdstep<float3>);
# endif
cuda::function<void (float4*, unsigned int)> lattice(mdsim::lattice<float3>);
cuda::function<void (float4*, unsigned int)> lattice_simple(mdsim::lattice_simple<float3>);
cuda::function<void (float4*, float, ushort3*)> boltzmann(mdsim::boltzmann);
#else /* DIM_3D */
cuda::function<void (float2*, float2*, float2*, float2 const*)> inteq(mdsim::inteq<float2>);
# ifdef USE_CELL
cuda::function<void (float2 const*, float2*, float2*, int const*, float*, float*)> mdstep(mdsim::mdstep<float2>);
cuda::function<void (float2 const*, int*)> assign_cells(mdsim::assign_cells<CELL_SIZE>);
cuda::function<void (int const*, int*)> update_neighbours(mdsim::update_neighbours<CELL_SIZE, float2>);
cuda::texture<float2> r(mdsim::t_r);
cuda::function<void (float2 const*, unsigned int*)> sfc_hilbert_encode(mdsim::sfc_hilbert_encode<float2>);
# else
cuda::function<void (float2*, float2*, float2*, float*, float*)> mdstep(mdsim::mdstep<float2>);
# endif
cuda::function<void (float2*, unsigned int)> lattice(mdsim::lattice<float2>);
cuda::function<void (float2*, unsigned int)> lattice_simple(mdsim::lattice_simple<float2>);
cuda::function<void (float2*, float, ushort3*)> boltzmann(mdsim::boltzmann);
#endif /* DIM_3D */

cuda::symbol<unsigned int> npart(mdsim::npart);
cuda::symbol<float> box(mdsim::box);
cuda::symbol<float> timestep(mdsim::timestep);
cuda::symbol<float> r_cut(mdsim::r_cut);
cuda::symbol<float> rr_cut(mdsim::rr_cut);
cuda::symbol<float> en_cut(mdsim::en_cut);
#ifdef USE_CELL
cuda::function<void (int const*, int*)> update_cells(mdsim::update_cells<CELL_SIZE>);
cuda::symbol<unsigned int> ncell(mdsim::ncell);
cuda::symbol<unsigned int> nnbl(mdsim::nnbl);
cuda::symbol<float> r_cell(mdsim::r_cell);
cuda::symbol<float> rr_cell(mdsim::rr_cell);
cuda::symbol<unsigned int> sfc_level(mdsim::sfc_level);
#endif

#ifdef USE_SMOOTH_POTENTIAL
cuda::symbol<float> rri_smooth(mdsim::rri_smooth);
cuda::function <void (float3*, const float2)> sample_smooth_function(sample_smooth_function);
#endif

cuda::symbol<uint3> a(::rand48::a);
cuda::symbol<uint3> c(::rand48::c);

}}} // namespace mdsim::gpu::ljfluid
