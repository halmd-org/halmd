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

#include <float.h>
#define LJFLUID_NAMESPACE cu_ljfluid_cell
#include <ljgpu/ljfluid/gpu/base.cuh>
#include <ljgpu/ljfluid/gpu/ljfluid_cell.hpp>

namespace ljgpu { namespace gpu { namespace LJFLUID_NAMESPACE
{

enum {
    /** fixed number of placeholders per cell */
    CELL_SIZE = ljfluid_base<ljfluid_impl_gpu_cell>::CELL_SIZE,
    /** virtual particle tag */
    VIRTUAL_PARTICLE = ljfluid_base<ljfluid_impl_gpu_cell>::VIRTUAL_PARTICLE,
};

/** number of cells per dimension */
static __constant__ uint ncell;

/**
 * determine cell index for a particle
 */
__device__ uint compute_cell(float3 r)
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
    r = (__saturatef(r / box) * (1.f - FLT_EPSILON)) * ncell;
    return (uint(r.z) * ncell + uint(r.y)) * ncell + uint(r.x);
}

__device__ uint compute_cell(float2 r)
{
    r = (__saturatef(r / box) * (1.f - FLT_EPSILON)) * ncell;
    return uint(r.y) * ncell + uint(r.x);
}

/**
 * compute neighbour cell
 */
__device__ uint compute_neighbour_cell(int3 const& offset)
{
    // cell belonging to this execution block
    int3 cell = make_int3(
		       blockIdx.x % ncell,
		       (blockIdx.x / ncell) % ncell,
		       blockIdx.x / ncell / ncell
		       );
    // neighbour cell of this cell
    int3 neighbour = make_int3(
			    (cell.x + ncell + offset.x) % ncell,
			    (cell.y + ncell + offset.y) % ncell,
			    (cell.z + ncell + offset.z) % ncell
			    );

    return (neighbour.z * ncell + neighbour.y) * ncell + neighbour.x;
}

__device__ uint compute_neighbour_cell(int2 const& offset)
{
    // cell belonging to this execution block
    int2 cell = make_int2(
		       blockIdx.x % ncell,
		       blockIdx.x / ncell
		       );
    // neighbour cell of this cell
    int2 neighbour = make_int2(
			    (cell.x + ncell + offset.x) % ncell,
			    (cell.y + ncell + offset.y) % ncell
			    );

    return neighbour.y * ncell + neighbour.x;
}

/**
 * compute forces with particles in a neighbour cell
 */
template <uint block_size, bool same_cell, typename T, typename TT, typename U, typename I>
__device__ void compute_cell_forces(U const* g_r, int const* g_n, I const& offset, T const& r, int const& n, TT& f, float& en, float& virial)
{
    __shared__ T s_r[block_size];
    __shared__ int s_n[block_size];

    // compute cell index
    uint cell = compute_neighbour_cell(offset);

    // load particles coordinates for cell
    s_r[threadIdx.x] = unpack(g_r[cell * block_size + threadIdx.x]);
    s_n[threadIdx.x] = g_n[cell * block_size + threadIdx.x];
    __syncthreads();

    if (n != VIRTUAL_PARTICLE) {
	for (uint i = 0; i < block_size; ++i) {
	    // skip placeholder particles
	    if (s_n[i] == VIRTUAL_PARTICLE)
		break;
	    // skip same particle
	    if (same_cell && threadIdx.x == i)
		continue;

	    compute_force(r, s_r[i], f, en, virial);
	}
    }
    __syncthreads();
}

/**
 * 3-dimensional MD simulation step
 */
template <uint block_size>
__global__ void mdstep(float4 const* g_r, float4* g_v, float4* g_f, int const* g_tag, float* g_en, float* g_virial)
{
    // load particle associated with this thread
    float3 r = unpack(g_r[GTID]);
    float3 v = unpack(g_v[GTID]);
    int tag = g_tag[GTID];

    // potential energy contribution
    float en = 0;
    // virial equation sum contribution
    float virial = 0;

#ifdef USE_DSFUN
    dfloat3 f(make_float3(0, 0, 0));
#else
    float3 f(make_float3(0, 0, 0));
#endif

#ifdef USE_CELL_SUMMATION_ORDER
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

    // sum forces over this cell
    compute_cell_forces<block_size, true>(g_r, g_tag, make_int3( 0,  0,  0), r, tag, f, en, virial);
    // sum forces over 26 neighbour cells, grouped into 13 pairs of mutually opposite cells
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3(-1, -1, -1), r, tag, f, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3(+1, +1, +1), r, tag, f, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3(-1, -1, +1), r, tag, f, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3(+1, +1, -1), r, tag, f, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3(-1, +1, +1), r, tag, f, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3(+1, -1, -1), r, tag, f, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3(+1, -1, +1), r, tag, f, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3(-1, +1, -1), r, tag, f, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3(-1, -1,  0), r, tag, f, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3(+1, +1,  0), r, tag, f, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3(-1, +1,  0), r, tag, f, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3(+1, -1,  0), r, tag, f, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3(-1,  0, -1), r, tag, f, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3(+1,  0, +1), r, tag, f, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3(-1,  0, +1), r, tag, f, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3(+1,  0, -1), r, tag, f, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3( 0, -1, -1), r, tag, f, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3( 0, +1, +1), r, tag, f, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3( 0, -1, +1), r, tag, f, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3( 0, +1, -1), r, tag, f, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3(-1,  0,  0), r, tag, f, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3(+1,  0,  0), r, tag, f, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3( 0, -1,  0), r, tag, f, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3( 0, +1,  0), r, tag, f, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3( 0,  0, -1), r, tag, f, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3( 0,  0, +1), r, tag, f, en, virial);

#else /* ! USE_CELL_SUMMATION_ORDER */
    // visit 26 neighbour cells
    compute_cell_forces<block_size, true>(g_r, g_tag, make_int3( 0,  0,  0), r, tag, f, en, virial);
    for (int x = -1; x <= 1; ++x)
	for (int y = -1; y <= 1; ++y)
	    for (int z = -1; z <= 1; ++z)
		if (x != 0 || y != 0 || z != 0)
		    compute_cell_forces<block_size, false>(g_r, g_tag, make_int3(x,  y,  z), r, tag, f, en, virial);
#endif /* USE_CELL_SUMMATION_ORDER */

    // second leapfrog step as part of integration of equations of motion
#ifdef USE_DSFUN
    leapfrog_full_step(v, f.f0);
#else
    leapfrog_full_step(v, f);
#endif

    // store particle associated with this thread
    g_v[GTID] = pack(v);
#ifdef USE_DSFUN
    g_f[GTID] = pack(f.f0);
#else
    g_f[GTID] = pack(f);
#endif
    g_en[GTID] = en;
    g_virial[GTID] = virial;
}

/**
 * 2-dimensional MD simulation step
 */
template <uint block_size>
__global__ void mdstep(float2 const* g_r, float2* g_v, float2* g_f, int const* g_tag, float* g_en, float* g_virial)
{
    // load particle associated with this thread
    float2 r = unpack(g_r[GTID]);
    float2 v = unpack(g_v[GTID]);
    int tag = g_tag[GTID];

    // potential energy contribution
    float en = 0;
    // virial equation sum contribution
    float virial = 0;

#ifdef USE_DSFUN
    dfloat2 f(make_float2(0, 0));
#else
    float2 f(make_float2(0, 0));
#endif

#ifdef USE_CELL_SUMMATION_ORDER
    // sum forces over this cell
    compute_cell_forces<block_size, true>(g_r, g_tag, make_int2( 0,  0), r, tag, f, en, virial);
    // sum forces over 8 neighbour cells, grouped into 4 pairs of mutually opposite cells
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int2(-1, -1), r, tag, f, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int2(+1, +1), r, tag, f, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int2(-1, +1), r, tag, f, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int2(+1, -1), r, tag, f, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int2(-1,  0), r, tag, f, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int2(+1,  0), r, tag, f, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int2( 0, -1), r, tag, f, en, virial);
    compute_cell_forces<block_size, false>(g_r, g_tag, make_int2( 0, +1), r, tag, f, en, virial);
#else
    compute_cell_forces<block_size, true>(g_r, g_tag, make_int2( 0,  0), r, tag, f, en, virial);
    // visit 8 neighbour cells
    for (int x = -1; x <= 1; ++x)
	for (int y = -1; y <= 1; ++y)
	    if (x != 0 || y != 0)
		compute_cell_forces<block_size, false>(g_r, g_tag, make_int2(x, y), r, tag, f, en, virial);
#endif

    // second leapfrog step as part of integration of equations of motion
#ifdef USE_DSFUN
    leapfrog_full_step(v, f.f0);
#else
    leapfrog_full_step(v, f);
#endif

    // store particle associated with this thread
    g_v[GTID] = pack(v);
#ifdef USE_DSFUN
    g_f[GTID] = pack(f.f0);
#else
    g_f[GTID] = pack(f);
#endif
    g_en[GTID] = en;
    g_virial[GTID] = virial;
}

/**
 * assign particles to cells
 */
template <uint cell_size, typename T, typename U>
__global__ void assign_cells(U const* g_part, U* g_r, int* g_tag)
{
    __shared__ T s_block[cell_size];
    __shared__ T s_r[cell_size];
    __shared__ int s_icell[cell_size];
    __shared__ int s_tag[cell_size];
    // number of particles in cell
    uint n = 0;

    // mark all particles in cell as virtual particles
    s_tag[threadIdx.x] = VIRTUAL_PARTICLE;

    __syncthreads();

    for (uint i = 0; i < npart; i += cell_size) {
	// load block of particles from global device memory
	T r = unpack(g_part[i + threadIdx.x]);
	s_block[threadIdx.x] = r;
	s_icell[threadIdx.x] = compute_cell(r);
	__syncthreads();

	if (threadIdx.x == 0) {
	    for (uint j = 0; j < cell_size && (i + j) < npart; j++) {
		if (s_icell[j] == blockIdx.x) {
		    // store particle in cell
		    s_r[n] = s_block[j];
		    // store particle number
		    s_tag[n] = i + j;
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
template <uint cell_size, typename T, typename U, typename I>
__device__ void examine_cell(I const& offset, U const* g_ir, U const* g_iR, U const* g_iv, int const* g_itag, T* s_or, T* s_oR, T* s_ov, int* s_otag, uint& npart)
{
    __shared__ T s_ir[cell_size];
    __shared__ T s_iR[cell_size];
    __shared__ T s_iv[cell_size];
    __shared__ int s_itag[cell_size];
    __shared__ uint s_cell[cell_size];

    // compute cell index
    uint cell = compute_neighbour_cell(offset);

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
	for (uint j = 0; j < cell_size; j++) {
	    // skip virtual particles
	    if (s_itag[j] == VIRTUAL_PARTICLE)
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

template <uint cell_size>
__device__ void examine_cells(float4 const* g_ir, float4 const* g_iR, float4 const* g_iv, int const* g_itag, float3* s_or, float3* s_oR, float3* s_ov, int* s_otag, uint& n)
{
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
}

template <uint cell_size>
__device__ void examine_cells(float2 const* g_ir, float2 const* g_iR, float2 const* g_iv, int const* g_itag, float2* s_or, float2* s_oR, float2* s_ov, int* s_otag, uint& n)
{
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
}

/**
 * update cells
 */
template <uint cell_size, typename T, typename U>
__global__ void update_cells(U const* g_ir, U const* g_iR, U const* g_iv, int const* g_itag, U* g_or, U* g_oR, U* g_ov, int* g_otag)
{
    __shared__ T s_or[cell_size];
    __shared__ T s_oR[cell_size];
    __shared__ T s_ov[cell_size];
    __shared__ int s_otag[cell_size];
    // number of particles in cell
    uint n = 0;

    // mark all particles in cell as virtual particles
    s_otag[threadIdx.x] = VIRTUAL_PARTICLE;
    __syncthreads();

    examine_cells<cell_size>(g_ir, g_iR, g_iv, g_itag, s_or, s_oR, s_ov, s_otag, n);

    // store cell in global device memory
    g_or[blockIdx.x * cell_size + threadIdx.x] = pack(s_or[threadIdx.x]);
    g_oR[blockIdx.x * cell_size + threadIdx.x] = pack(s_oR[threadIdx.x]);
    g_ov[blockIdx.x * cell_size + threadIdx.x] = pack(s_ov[threadIdx.x]);
    g_otag[blockIdx.x * cell_size + threadIdx.x] = s_otag[threadIdx.x];
}

} // namespace ljgpu::gpu::LJFLUID_NAMESPACE

/**
 * device constant wrappers
 */
cuda::symbol<uint> ljfluid_base<ljfluid_impl_gpu_cell>::npart(LJFLUID_NAMESPACE::npart);
cuda::symbol<float> ljfluid_base<ljfluid_impl_gpu_cell>::box(LJFLUID_NAMESPACE::box);
cuda::symbol<float> ljfluid_base<ljfluid_impl_gpu_cell>::timestep(LJFLUID_NAMESPACE::timestep);
cuda::symbol<float> ljfluid_base<ljfluid_impl_gpu_cell>::r_cut(LJFLUID_NAMESPACE::r_cut);
cuda::symbol<float> ljfluid_base<ljfluid_impl_gpu_cell>::rr_cut(LJFLUID_NAMESPACE::rr_cut);
cuda::symbol<float> ljfluid_base<ljfluid_impl_gpu_cell>::en_cut(LJFLUID_NAMESPACE::en_cut);
cuda::symbol<float> ljfluid_base<ljfluid_impl_gpu_cell>::rri_smooth(LJFLUID_NAMESPACE::rri_smooth);
cuda::symbol<uint> ljfluid_base<ljfluid_impl_gpu_cell>::ncell(LJFLUID_NAMESPACE::ncell);

/**
 * device function wrappers
 */
cuda::function<void (float3*, const float2)>
    ljfluid_base<ljfluid_impl_gpu_cell>::sample_smooth_function(LJFLUID_NAMESPACE::sample_smooth_function);

cuda::function<void (float4*, float4*, float4*, float4 const*)>
    ljfluid<ljfluid_impl_gpu_cell<3> >::inteq(LJFLUID_NAMESPACE::inteq<float3>);
cuda::function<void (float4 const*, float4*, float4*, int const*, float*, float*)>
    ljfluid<ljfluid_impl_gpu_cell<3> >::mdstep(LJFLUID_NAMESPACE::mdstep<CELL_SIZE>);
cuda::function<void (float4 const*, float4*, int*)>
    ljfluid<ljfluid_impl_gpu_cell<3> >::assign_cells(LJFLUID_NAMESPACE::assign_cells<CELL_SIZE, float3>);
cuda::function<void (float4 const*, float4 const*, float4 const*, int const*, float4*, float4*, float4*, int*)>
    ljfluid<ljfluid_impl_gpu_cell<3> >::update_cells(LJFLUID_NAMESPACE::update_cells<CELL_SIZE, float3>);

cuda::function<void (float2*, float2*, float2*, float2 const*)>
    ljfluid<ljfluid_impl_gpu_cell<2> >::inteq(LJFLUID_NAMESPACE::inteq<float2>);
cuda::function<void (float2 const*, float2*, float2*, int const*, float*, float*)>
    ljfluid<ljfluid_impl_gpu_cell<2> >::mdstep(LJFLUID_NAMESPACE::mdstep<CELL_SIZE>);
cuda::function<void (float2 const*, float2*, int*)>
    ljfluid<ljfluid_impl_gpu_cell<2> >::assign_cells(LJFLUID_NAMESPACE::assign_cells<CELL_SIZE, float2>);
cuda::function<void (float2 const*, float2 const*, float2 const*, int const*, float2*, float2*, float2*, int*)>
    ljfluid<ljfluid_impl_gpu_cell<2> >::update_cells(LJFLUID_NAMESPACE::update_cells<CELL_SIZE, float2>);

}} // namespace ljgpu::gpu
