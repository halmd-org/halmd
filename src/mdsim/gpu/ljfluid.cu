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
static __constant__ unsigned int nbl_size;
/** neighbour list stride */
static __constant__ unsigned int nbl_stride;
/** potential cutoff distance with cell skin */
static __constant__ float r_cell;
/** squared potential cutoff distance with cell skin */
static __constant__ float rr_cell;
# ifdef DIM_3D
/** texture reference to periodic particle positions */
struct texture<float4, 1, cudaReadModeElementType> t_r;
/** texture reference to extended particle positions */
struct texture<float4, 1, cudaReadModeElementType> t_R;
/** texture reference to particle velocities */
struct texture<float4, 1, cudaReadModeElementType> t_v;
# else
/** texture reference to periodic particle positions */
struct texture<float2, 1, cudaReadModeElementType> t_r;
/** texture reference to extended particle positions */
struct texture<float2, 1, cudaReadModeElementType> t_R;
/** texture reference to particle velocities */
struct texture<float2, 1, cudaReadModeElementType> t_v;
# endif
/** texture reference to particle tags */
struct texture<int, 1, cudaReadModeElementType> t_tag;
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

    for (unsigned int i = 0; i < nbl_size; ++i) {
	// coalesced read from neighbour list
	const int n = g_nbl[i * nbl_stride + GTID];
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

/**
 * blockwise maximum velocity magnitude
 */
template <typename T, typename U>
__global__ void maximum_velocity(U const* g_v, float* g_vmax)
{
    extern __shared__ float s_vv[];

    // load first particle from global device memory
    T v = unpack(g_v[GTID]);
    float vv = v * v;
    // load further particles
    for (unsigned int i = GTDIM + GTID; i < npart; i += GTDIM) {
	v = unpack(g_v[i]);
	vv = fmaxf(vv, v * v);
    }
    // maximum velocity for this thread
    s_vv[TID] = vv;
    __syncthreads();

    // compute maximum velocity for all threads in block
    if (TID < 256) {
	vv = fmaxf(vv, s_vv[TID + 256]);
	s_vv[TID] = vv;
    }
    __syncthreads();
    if (TID < 128) {
	vv = fmaxf(vv, s_vv[TID + 128]);
	s_vv[TID] = vv;
    }
    __syncthreads();
    if (TID < 64) {
	vv = fmaxf(vv, s_vv[TID + 64]);
	s_vv[TID] = vv;
    }
    __syncthreads();
    if (TID < 32) {
	vv = fmaxf(vv, s_vv[TID + 32]);
	s_vv[TID] = vv;
    }
    // no further syncs needed within execution warp of 32 threads
    if (TID < 16) {
	vv = fmaxf(vv, s_vv[TID + 16]);
	s_vv[TID] = vv;
    }
    if (TID < 8) {
	vv = fmaxf(vv, s_vv[TID + 8]);
	s_vv[TID] = vv;
    }
    if (TID < 4) {
	vv = fmaxf(vv, s_vv[TID + 4]);
	s_vv[TID] = vv;
    }
    if (TID < 2) {
	vv = fmaxf(vv, s_vv[TID + 2]);
	s_vv[TID] = vv;
    }
    if (TID < 1) {
	vv = fmaxf(vv, s_vv[TID + 1]);
	// store maximum block velocity in global memory
	g_vmax[blockIdx.x] = sqrtf(vv);
    }
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
 * blockwise potential energy sum
 */
__global__ void potential_energy_sum(float const* g_en, float2* g_en_sum)
{
    // single-double floating point arithmetic
    extern __shared__ dfloat s_en[];

    // load first particle from global device memory
    dfloat en = g_en[GTID];
    // load further particles
    for (unsigned int i = GTDIM + GTID; i < npart; i += GTDIM) {
	en = en + g_en[i];
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
 * initialise particle tags
 */
__global__ void init_tags(int* g_tag)
{
    int tag = VIRTUAL_PARTICLE;
    if (GTID < npart) {
	tag = GTID;
    }
    g_tag[GTID] = tag;
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
	    if (rr <= rr_cell && count < nbl_size) {
		// scattered write to neighbour list
		g_nbl[count * nbl_stride + n] = m;
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
}

/**
 * compute cell indices for given particle positions
 */
template <typename T, typename U>
__global__ void compute_cell(U const* g_part, uint* g_cell)
{
    const T r = unpack(g_part[GTID]);

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
    g_cell[GTID] = uint(cell.x) + ncell * (uint(cell.y) + ncell * uint(cell.z));
#else
    g_cell[GTID] = uint(cell.x) + ncell * uint(cell.y);
#endif
}

/**
 * compute global cell offsets in particle list
 */
__global__ void find_cell_offset(uint* g_cell, int* g_cell_offset)
{
    const uint j = g_cell[GTID];
    const uint k = (GTID > 0 && GTID < npart) ? g_cell[GTID - 1] : j;

    if (GTID == 0 || k < j) {
	// particle marks the start of a cell
	g_cell_offset[j] = GTID;
    }
}

/**
 * assign particles to cells
 */
template <uint cell_size>
__global__ void assign_cells(uint const* g_cell, int const* g_cell_offset, int const* g_itag, int* g_otag)
{
    __shared__ int s_offset[1];

    if (threadIdx.x == 0) {
	// global offset of this cell in particle list
	s_offset[0] = g_cell_offset[blockIdx.x];
    }
    __syncthreads();

    const int offset = s_offset[0];
    // mark as virtual particle
    int tag = -1;

    if (offset >= 0) {
	const int n = offset + threadIdx.x;
	const uint cell = g_cell[n];
	const uint itag = g_itag[n];

	// assign particle to cell
	if (n < npart && cell == blockIdx.x) {
	    tag = itag;
	}
    }

    // store cell in global device memory
    g_otag[blockIdx.x * cell_size + threadIdx.x] = tag;
}

/**
 * generate ascending index sequence
 */
__global__ void gen_index(int* g_idx)
{
    g_idx[GTID] = (GTID < npart) ? GTID : 0;
}

/**
 * order particles after given permutation
 */
template <typename U>
__global__ void order_particles(const int* g_idx, U* g_or, U* g_oR, U* g_ov, int* g_otag)
{
    // permutation index
    const uint j = g_idx[GTID];
    // permute particle phase space coordinates
    g_or[GTID] = tex1Dfetch(t_r, j);
    g_oR[GTID] = tex1Dfetch(t_R, j);
    g_ov[GTID] = tex1Dfetch(t_v, j);
    // permute particle tracking number
    g_otag[GTID] = tex1Dfetch(t_tag, j);
}

#endif  /* USE_CELL */

} // namespace mdsim


namespace mdsim { namespace gpu { namespace ljfluid
{

#ifdef DIM_3D
cuda::function<void (float4*, float4*, float4*, float4 const*)> inteq(mdsim::inteq<float3>);
# ifdef USE_CELL
cuda::function<void (float4 const*, float4*, float4*, int const*, float*, float*)> mdstep(mdsim::mdstep<float3>);
cuda::function<void (float4 const*, float*)> maximum_velocity(mdsim::maximum_velocity<float3>);
cuda::function<void (int const*, int*)> update_neighbours(mdsim::update_neighbours<CELL_SIZE, float3>);
cuda::function<void (float4 const*, uint*)> compute_cell(mdsim::compute_cell<float3>);
cuda::function<void (const int*, float4*, float4*, float4*, int*)> order_particles(mdsim::order_particles);
cuda::texture<float4> r(mdsim::t_r);
cuda::texture<float4> R(mdsim::t_R);
cuda::texture<float4> v(mdsim::t_v);
# else
cuda::function<void (float4*, float4*, float4*, float*, float*)> mdstep(mdsim::mdstep<float3>);
# endif
cuda::function<void (float4*, unsigned int)> lattice(mdsim::lattice<float3>);
cuda::function<void (float4*, unsigned int)> lattice_simple(mdsim::lattice_simple<float3>);
#else /* DIM_3D */
cuda::function<void (float2*, float2*, float2*, float2 const*)> inteq(mdsim::inteq<float2>);
# ifdef USE_CELL
cuda::function<void (float2 const*, float2*, float2*, int const*, float*, float*)> mdstep(mdsim::mdstep<float2>);
cuda::function<void (float2 const*, float*)> maximum_velocity(mdsim::maximum_velocity<float2>);
cuda::function<void (int const*, int*)> update_neighbours(mdsim::update_neighbours<CELL_SIZE, float2>);
cuda::function<void (float2 const*, uint*)> compute_cell(mdsim::compute_cell<float2>);
cuda::function<void (const int*, float2*, float2*, float2*, int*)> order_particles(mdsim::order_particles);
cuda::texture<float2> r(mdsim::t_r);
cuda::texture<float2> R(mdsim::t_R);
cuda::texture<float2> v(mdsim::t_v);
# else
cuda::function<void (float2*, float2*, float2*, float*, float*)> mdstep(mdsim::mdstep<float2>);
# endif
cuda::function<void (float2*, unsigned int)> lattice(mdsim::lattice<float2>);
cuda::function<void (float2*, unsigned int)> lattice_simple(mdsim::lattice_simple<float2>);
#endif /* DIM_3D */

cuda::symbol<unsigned int> npart(mdsim::npart);
cuda::symbol<float> box(mdsim::box);
cuda::symbol<float> timestep(mdsim::timestep);
cuda::symbol<float> r_cut(mdsim::r_cut);
cuda::symbol<float> rr_cut(mdsim::rr_cut);
cuda::symbol<float> en_cut(mdsim::en_cut);
cuda::function<void (int*)> init_tags(mdsim::init_tags);
cuda::function<void (float const* g_en, float2* g_en_sum)> potential_energy_sum(mdsim::potential_energy_sum);

#ifdef USE_CELL
cuda::symbol<unsigned int> ncell(mdsim::ncell);
cuda::symbol<unsigned int> nbl_size(mdsim::nbl_size);
cuda::symbol<unsigned int> nbl_stride(mdsim::nbl_stride);
cuda::symbol<float> r_cell(mdsim::r_cell);
cuda::symbol<float> rr_cell(mdsim::rr_cell);
cuda::texture<int> tag(mdsim::t_tag);
cuda::function<void (uint const*, int const*, int const*, int*)> assign_cells(mdsim::assign_cells<CELL_SIZE>);
cuda::function<void (uint*, int*)> find_cell_offset(mdsim::find_cell_offset);
cuda::function<void (int*)> gen_index(mdsim::gen_index);
#endif

#ifdef USE_SMOOTH_POTENTIAL
cuda::symbol<float> rri_smooth(mdsim::rri_smooth);
cuda::function <void (float3*, const float2)> sample_smooth_function(sample_smooth_function);
#endif

}}} // namespace mdsim::gpu::ljfluid
