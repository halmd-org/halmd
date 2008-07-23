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
#ifdef DIM_3D
struct texture<float4, 1, cudaReadModeElementType> t_r;
#else
struct texture<float2, 1, cudaReadModeElementType> t_r;
#endif
struct texture<int, 1, cudaReadModeElementType> t_cell;
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
#ifdef DIM_3D
__device__ int3 compute_cell(T const& r)
#else
__device__ int2 compute_cell(T const& r)
#endif
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
    return make_int3(cell.x, cell.y, cell.z);
#else
    return make_int2(cell.x, cell.y);
#endif
}

template <typename T>
__device__ void compute_cell(T const& r, unsigned int& i)
{
#ifdef DIM_3D
    const int3 cell = compute_cell(r);
    i = cell.x + ncell * (cell.y + ncell * cell.z);
#else
    const int2 cell = compute_cell(r);
    i = cell.x + ncell * cell.y;
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
 * compute neighbour cell
 */
template <typename I>
__device__ unsigned int compute_neighbour_cell(I const& cell, I const& offset)
{
#ifdef DIM_3D
    I neighbour = make_int3((cell.x + ncell + offset.x) % ncell, (cell.y + ncell + offset.y) % ncell, (cell.z + ncell + offset.z) % ncell);
    return (neighbour.z * ncell + neighbour.y) * ncell + neighbour.x;
#else
    I neighbour = make_int2((cell.x + ncell + offset.x) % ncell, (cell.y + ncell + offset.y) % ncell);
    return neighbour.y * ncell + neighbour.x;
#endif
}

/**
 * compute forces with particles in a neighbour cell
 */
template <unsigned int block_size, typename T, typename TT, typename U, typename I>
__device__ void compute_cell_forces(U const* g_r, I const& offset, T const& r, I const& cell, TT& f, float& en, float& virial)
{
    // neighbour cell index
    const unsigned int neighbour = compute_neighbour_cell(cell, offset);

    for (unsigned int i = 0; i < block_size; ++i) {
	// load particle number from cell placeholder
	const int n = tex1Dfetch(t_cell, neighbour * block_size + i);
	// skip placeholder particles
	if (!IS_REAL_PARTICLE(n))
	    break;
	// skip same particle
	if (GTID == n)
	    continue;

	compute_force(r, unpack(tex1Dfetch(t_r, n)), f, en, virial);
    }
}

/**
 * n-dimensional MD simulation step
 */
template <unsigned int block_size, typename T, typename U>
__global__ void mdstep(U const* g_r, U* g_v, U* g_f, float* g_en, float* g_virial)
{
    // load particle associated with this thread
    const T r = unpack(g_r[GTID]);
    T v = unpack(g_v[GTID]);

#ifdef DIM_3D
    // cell to which particle belongs
    const int3 cell = compute_cell(r);
#else
    const int2 cell = compute_cell(r);
#endif

    // potential energy contribution
    float en = 0;
    // virial equation sum contribution
    float virial = 0;

#ifdef DIM_3D
    dfloat3 f(make_float3(0, 0, 0));
#else
    dfloat2 f(make_float2(0, 0));
#endif

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
    // sum forces over this cell
    compute_cell_forces<block_size>(g_r, make_int3( 0,  0,  0), r, cell, f, en, virial);
    // sum forces over 26 neighbour cells, grouped into 13 pairs of mutually opposite cells
    compute_cell_forces<block_size>(g_r, make_int3(-1, -1, -1), r, cell, f, en, virial);
    compute_cell_forces<block_size>(g_r, make_int3(+1, +1, +1), r, cell, f, en, virial);
    compute_cell_forces<block_size>(g_r, make_int3(-1, -1, +1), r, cell, f, en, virial);
    compute_cell_forces<block_size>(g_r, make_int3(+1, +1, -1), r, cell, f, en, virial);
    compute_cell_forces<block_size>(g_r, make_int3(-1, +1, +1), r, cell, f, en, virial);
    compute_cell_forces<block_size>(g_r, make_int3(+1, -1, -1), r, cell, f, en, virial);
    compute_cell_forces<block_size>(g_r, make_int3(+1, -1, +1), r, cell, f, en, virial);
    compute_cell_forces<block_size>(g_r, make_int3(-1, +1, -1), r, cell, f, en, virial);
    compute_cell_forces<block_size>(g_r, make_int3(-1, -1,  0), r, cell, f, en, virial);
    compute_cell_forces<block_size>(g_r, make_int3(+1, +1,  0), r, cell, f, en, virial);
    compute_cell_forces<block_size>(g_r, make_int3(-1, +1,  0), r, cell, f, en, virial);
    compute_cell_forces<block_size>(g_r, make_int3(+1, -1,  0), r, cell, f, en, virial);
    compute_cell_forces<block_size>(g_r, make_int3(-1,  0, -1), r, cell, f, en, virial);
    compute_cell_forces<block_size>(g_r, make_int3(+1,  0, +1), r, cell, f, en, virial);
    compute_cell_forces<block_size>(g_r, make_int3(-1,  0, +1), r, cell, f, en, virial);
    compute_cell_forces<block_size>(g_r, make_int3(+1,  0, -1), r, cell, f, en, virial);
    compute_cell_forces<block_size>(g_r, make_int3( 0, -1, -1), r, cell, f, en, virial);
    compute_cell_forces<block_size>(g_r, make_int3( 0, +1, +1), r, cell, f, en, virial);
    compute_cell_forces<block_size>(g_r, make_int3( 0, -1, +1), r, cell, f, en, virial);
    compute_cell_forces<block_size>(g_r, make_int3( 0, +1, -1), r, cell, f, en, virial);
    compute_cell_forces<block_size>(g_r, make_int3(-1,  0,  0), r, cell, f, en, virial);
    compute_cell_forces<block_size>(g_r, make_int3(+1,  0,  0), r, cell, f, en, virial);
    compute_cell_forces<block_size>(g_r, make_int3( 0, -1,  0), r, cell, f, en, virial);
    compute_cell_forces<block_size>(g_r, make_int3( 0, +1,  0), r, cell, f, en, virial);
    compute_cell_forces<block_size>(g_r, make_int3( 0,  0, -1), r, cell, f, en, virial);
    compute_cell_forces<block_size>(g_r, make_int3( 0,  0, +1), r, cell, f, en, virial);
#else
    // sum forces over this cell
    compute_cell_forces<block_size>(g_r, make_int2( 0,  0), r, cell, f, en, virial);
    // sum forces over 8 neighbour cells, grouped into 4 pairs of mutually opposite cells
    compute_cell_forces<block_size>(g_r, make_int2(-1, -1), r, cell, f, en, virial);
    compute_cell_forces<block_size>(g_r, make_int2(+1, +1), r, cell, f, en, virial);
    compute_cell_forces<block_size>(g_r, make_int2(-1, +1), r, cell, f, en, virial);
    compute_cell_forces<block_size>(g_r, make_int2(+1, -1), r, cell, f, en, virial);
    compute_cell_forces<block_size>(g_r, make_int2(-1,  0), r, cell, f, en, virial);
    compute_cell_forces<block_size>(g_r, make_int2(+1,  0), r, cell, f, en, virial);
    compute_cell_forces<block_size>(g_r, make_int2( 0, -1), r, cell, f, en, virial);
    compute_cell_forces<block_size>(g_r, make_int2( 0, +1), r, cell, f, en, virial);
#endif

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
	compute_cell(unpack(g_part[i + threadIdx.x]), s_cell[threadIdx.x]);
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

#ifdef DIM_3D
    // compute cell index
    const unsigned int cell = compute_neighbour_cell(make_int3(blockIdx.x % ncell, (blockIdx.x / ncell) % ncell, blockIdx.x / ncell / ncell), offset);
#else
    const unsigned int cell = compute_neighbour_cell(make_int2(blockIdx.x % ncell, blockIdx.x / ncell), offset);
#endif

    // load particle numbers from global device memory
    int n = g_itag[cell * cell_size + threadIdx.x];
    s_itag[threadIdx.x] = n;

    if (IS_REAL_PARTICLE(n)) {
	// compute new cell
	compute_cell(unpack(tex1Dfetch(t_r, n)), s_cell[threadIdx.x]);
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
    examine_cell<cell_size>(make_int3( 0,  0,  0), g_itag, s_otag, n);
    // visit 26 neighbour cells
    examine_cell<cell_size>(make_int3(-1,  0,  0), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3(+1,  0,  0), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3( 0, -1,  0), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3( 0, +1,  0), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3(-1, -1,  0), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3(-1, +1,  0), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3(+1, -1,  0), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3(+1, +1,  0), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3( 0,  0, -1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3(-1,  0, -1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3(+1,  0, -1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3( 0, -1, -1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3( 0, +1, -1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3(-1, -1, -1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3(-1, +1, -1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3(+1, -1, -1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3(+1, +1, -1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3( 0,  0, +1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3(-1,  0, +1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3(+1,  0, +1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3( 0, -1, +1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3( 0, +1, +1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3(-1, -1, +1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3(-1, +1, +1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3(+1, -1, +1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int3(+1, +1, +1), g_itag, s_otag, n);
#else
    examine_cell<cell_size>(make_int2( 0,  0), g_itag, s_otag, n);
    // visit 8 neighbour cells
    examine_cell<cell_size>(make_int2(-1,  0), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int2(+1,  0), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int2( 0, -1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int2( 0, +1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int2(-1, -1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int2(-1, +1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int2(+1, -1), g_itag, s_otag, n);
    examine_cell<cell_size>(make_int2(+1, +1), g_itag, s_otag, n);
#endif

    // store cell in global device memory
    g_otag[blockIdx.x * cell_size + threadIdx.x] = s_otag[threadIdx.x];
}
#endif  /* USE_CELL */

} // namespace mdsim


namespace mdsim { namespace gpu { namespace ljfluid
{

#ifdef DIM_3D
cuda::function<void (float4*, float4*, float4*, float4 const*)> inteq(mdsim::inteq<float3>);
#ifdef USE_CELL
cuda::function<void (float4 const*, float4*, float4*, float*, float*)> mdstep(mdsim::mdstep<CELL_SIZE, float3>);
cuda::function<void (float4 const*, int*)> assign_cells(mdsim::assign_cells<CELL_SIZE>);
cuda::function<void (int const*, int*)> update_cells(mdsim::update_cells<CELL_SIZE>);
cuda::texture<float4> r(mdsim::t_r);
#else
cuda::function<void (float4*, float4*, float4*, float*, float*)> mdstep(mdsim::mdstep<float3>);
#endif
cuda::function<void (float4*, unsigned int)> lattice(mdsim::lattice<float3>);
cuda::function<void (float4*, float, ushort3*)> boltzmann(mdsim::boltzmann);
#else /* DIM_3D */
cuda::function<void (float2*, float2*, float2*, float2 const*)> inteq(mdsim::inteq<float2>);
#ifdef USE_CELL
cuda::function<void (float2 const*, float2*, float2*, float*, float*)> mdstep(mdsim::mdstep<CELL_SIZE, float2>);
cuda::function<void (float2 const*, int*)> assign_cells(mdsim::assign_cells<CELL_SIZE>);
cuda::function<void (int const*, int*)> update_cells(mdsim::update_cells<CELL_SIZE>);
cuda::texture<float2> r(mdsim::t_r);
#else
cuda::function<void (float2*, float2*, float2*, float*, float*)> mdstep(mdsim::mdstep<float2>);
#endif
cuda::function<void (float2*, unsigned int)> lattice(mdsim::lattice<float2>);
cuda::function<void (float2*, float, ushort3*)> boltzmann(mdsim::boltzmann);
#endif /* DIM_3D */

cuda::symbol<unsigned int> npart(mdsim::npart);
cuda::symbol<float> box(mdsim::box);
cuda::symbol<float> timestep(mdsim::timestep);
cuda::symbol<float> r_cut(mdsim::r_cut);
cuda::symbol<float> rr_cut(mdsim::rr_cut);
cuda::symbol<float> en_cut(mdsim::en_cut);
#ifdef USE_CELL
cuda::symbol<unsigned int> ncell(mdsim::ncell);
cuda::texture<int> cell(mdsim::t_cell);
#endif

#ifdef USE_SMOOTH_POTENTIAL
cuda::symbol<float> rri_smooth(mdsim::rri_smooth);
cuda::function <void (float3*, const float2)> sample_smooth_function(sample_smooth_function);
#endif

cuda::symbol<uint3> a(::rand48::a);
cuda::symbol<uint3> c(::rand48::c);

}}} // namespace mdsim::gpu::ljfluid
