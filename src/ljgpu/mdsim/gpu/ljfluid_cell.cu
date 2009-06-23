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
#include <ljgpu/mdsim/gpu/base.cuh>
#include <ljgpu/mdsim/gpu/ljfluid_cell.hpp>
using namespace ljgpu::gpu;

namespace ljgpu { namespace cu { namespace ljfluid
{

/** fixed number of placeholders per cell */
enum { CELL_SIZE = ljfluid_base<ljfluid_impl_gpu_cell>::CELL_SIZE };

/** number of cells per dimension */
__constant__ uint ncell;

/**
 * determine cell index for a particle
 */
__device__ uint compute_cell(vector<float, 3> r)
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

__device__ uint compute_cell(vector<float, 2> r)
{
    r = (__saturatef(r / box) * (1.f - FLT_EPSILON)) * ncell;
    return uint(r.y) * ncell + uint(r.x);
}

/**
 */
__device__ uint compute_neighbour_cell(int3 const& offset)
{
    // cell belonging to this execution block
    int x = BID % ncell;
    int y = (BID / ncell) % ncell;
    int z = BID / ncell / ncell;
    // neighbour cell of this cell
    x = (x + ncell + offset.x) % ncell;
    y = (y + ncell + offset.y) % ncell;
    z = (z + ncell + offset.z) % ncell;

    return (z * ncell + y) * ncell + x;
}

__device__ uint compute_neighbour_cell(int2 const& offset)
{
    // cell belonging to this execution block
    int x = BID % ncell;
    int y = BID / ncell;
    // neighbour cell of this cell
    x = (x + ncell + offset.x) % ncell;
    y = (y + ncell + offset.y) % ncell;

    return y * ncell + x;
}

/**
 * compute forces with particles in a neighbour cell
 */
template <bool same_cell,
	  mixture_type mixture,
	  potential_type potential,
	  typename I, typename T, typename U, typename V>
__device__ void compute_cell_forces(float4 const* g_r, I const& offset,
				    T const& r, unsigned int const tag, U& f,
				    float& en, V& virial)
{
    __shared__ T s_r[CELL_SIZE];
    __shared__ unsigned int s_tag[CELL_SIZE];

    // shared memory barrier
    __syncthreads();

    // compute cell index
    uint cell = compute_neighbour_cell(offset);
    // load particles coordinates for cell
    s_r[threadIdx.x] = detach_particle_tag(g_r[cell * CELL_SIZE + threadIdx.x], s_tag[threadIdx.x]);
    __syncthreads();

    if (tag == VIRTUAL_PARTICLE) return;

    // particle type in binary mixture
    int const a = (tag >= mpart[0]);

    for (uint i = 0; i < CELL_SIZE; ++i) {
	// skip placeholder particles
	if (s_tag[i] == VIRTUAL_PARTICLE) break;
	// skip same particle
	if (same_cell && threadIdx.x == i) continue;

	// particle type in binary mixture
	int const b = (s_tag[i] >= mpart[0]);

	compute_force<mixture, potential>(r, s_r[i], f, en, virial, a + b);
    }
}

/**
 * 3-dimensional MD simulation step
 */
template <typename vector_type,
	  mixture_type mixture,
	  potential_type potential,
	  typename T>
__global__ void mdstep(float4 const* g_r, T* g_v, T* g_f, float* g_en, T* g_virial)
{
    enum { dimension = vector_type::static_size };

    // load particle associated with this thread
    unsigned int tag;
    vector_type r = detach_particle_tag(g_r[GTID], tag);
    vector_type v = g_v[GTID];
    // potential energy contribution
    float en = 0;
    // virial equation sum contribution
    vector<float, (dimension - 1) * dimension / 2 + 1> virial = 0;
    // force sum
#ifdef USE_FORCE_DSFUN
    vector<dsfloat, dimension> f = 0;
#else
    vector_type f = 0;
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

    if (dimension == 3) {
	// sum forces over this cell
	compute_cell_forces<true, mixture, potential>(g_r, make_int3( 0,  0,  0), r, tag, f, en, virial);
	// sum forces over 26 neighbour cells, grouped into 13 pairs of mutually opposite cells
	compute_cell_forces<false, mixture, potential>(g_r, make_int3(-1, -1, -1), r, tag, f, en, virial);
	compute_cell_forces<false, mixture, potential>(g_r, make_int3(+1, +1, +1), r, tag, f, en, virial);
	compute_cell_forces<false, mixture, potential>(g_r, make_int3(-1, -1, +1), r, tag, f, en, virial);
	compute_cell_forces<false, mixture, potential>(g_r, make_int3(+1, +1, -1), r, tag, f, en, virial);
	compute_cell_forces<false, mixture, potential>(g_r, make_int3(-1, +1, +1), r, tag, f, en, virial);
	compute_cell_forces<false, mixture, potential>(g_r, make_int3(+1, -1, -1), r, tag, f, en, virial);
	compute_cell_forces<false, mixture, potential>(g_r, make_int3(+1, -1, +1), r, tag, f, en, virial);
	compute_cell_forces<false, mixture, potential>(g_r, make_int3(-1, +1, -1), r, tag, f, en, virial);
	compute_cell_forces<false, mixture, potential>(g_r, make_int3(-1, -1,  0), r, tag, f, en, virial);
	compute_cell_forces<false, mixture, potential>(g_r, make_int3(+1, +1,  0), r, tag, f, en, virial);
	compute_cell_forces<false, mixture, potential>(g_r, make_int3(-1, +1,  0), r, tag, f, en, virial);
	compute_cell_forces<false, mixture, potential>(g_r, make_int3(+1, -1,  0), r, tag, f, en, virial);
	compute_cell_forces<false, mixture, potential>(g_r, make_int3(-1,  0, -1), r, tag, f, en, virial);
	compute_cell_forces<false, mixture, potential>(g_r, make_int3(+1,  0, +1), r, tag, f, en, virial);
	compute_cell_forces<false, mixture, potential>(g_r, make_int3(-1,  0, +1), r, tag, f, en, virial);
	compute_cell_forces<false, mixture, potential>(g_r, make_int3(+1,  0, -1), r, tag, f, en, virial);
	compute_cell_forces<false, mixture, potential>(g_r, make_int3( 0, -1, -1), r, tag, f, en, virial);
	compute_cell_forces<false, mixture, potential>(g_r, make_int3( 0, +1, +1), r, tag, f, en, virial);
	compute_cell_forces<false, mixture, potential>(g_r, make_int3( 0, -1, +1), r, tag, f, en, virial);
	compute_cell_forces<false, mixture, potential>(g_r, make_int3( 0, +1, -1), r, tag, f, en, virial);
	compute_cell_forces<false, mixture, potential>(g_r, make_int3(-1,  0,  0), r, tag, f, en, virial);
	compute_cell_forces<false, mixture, potential>(g_r, make_int3(+1,  0,  0), r, tag, f, en, virial);
	compute_cell_forces<false, mixture, potential>(g_r, make_int3( 0, -1,  0), r, tag, f, en, virial);
	compute_cell_forces<false, mixture, potential>(g_r, make_int3( 0, +1,  0), r, tag, f, en, virial);
	compute_cell_forces<false, mixture, potential>(g_r, make_int3( 0,  0, -1), r, tag, f, en, virial);
	compute_cell_forces<false, mixture, potential>(g_r, make_int3( 0,  0, +1), r, tag, f, en, virial);
    }
    else {
	// sum forces over this cell
	compute_cell_forces<true, mixture, potential>(g_r, make_int2( 0,  0), r, tag, f, en, virial);
	// sum forces over 8 neighbour cells, grouped into 4 pairs of mutually opposite cells
	compute_cell_forces<false, mixture, potential>(g_r, make_int2(-1, -1), r, tag, f, en, virial);
	compute_cell_forces<false, mixture, potential>(g_r, make_int2(+1, +1), r, tag, f, en, virial);
	compute_cell_forces<false, mixture, potential>(g_r, make_int2(-1, +1), r, tag, f, en, virial);
	compute_cell_forces<false, mixture, potential>(g_r, make_int2(+1, -1), r, tag, f, en, virial);
	compute_cell_forces<false, mixture, potential>(g_r, make_int2(-1,  0), r, tag, f, en, virial);
	compute_cell_forces<false, mixture, potential>(g_r, make_int2(+1,  0), r, tag, f, en, virial);
	compute_cell_forces<false, mixture, potential>(g_r, make_int2( 0, -1), r, tag, f, en, virial);
	compute_cell_forces<false, mixture, potential>(g_r, make_int2( 0, +1), r, tag, f, en, virial);
    }

#else /* ! USE_CELL_SUMMATION_ORDER */
    if (dimension == 3) {
	compute_cell_forces<true, mixture, potential>(g_r, make_int3( 0,  0,  0), r, tag, f, en, virial);
	// visit 26 neighbour cells
	for (int x = -1; x <= 1; ++x)
	    for (int y = -1; y <= 1; ++y)
		for (int z = -1; z <= 1; ++z)
		    if (x != 0 || y != 0 || z != 0)
			compute_cell_forces<false, mixture, potential>(g_r, make_int3(x,  y,  z), r, tag, f, en, virial);
    }
    else {
	compute_cell_forces<true, mixture, potential>(g_r, make_int2( 0,  0), r, tag, f, en, virial);
	// visit 8 neighbour cells
	for (int x = -1; x <= 1; ++x)
	    for (int y = -1; y <= 1; ++y)
		if (x != 0 || y != 0)
		    compute_cell_forces<false, mixture, potential>(g_r, make_int2(x, y), r, tag, f, en, virial);
    }
#endif /* USE_CELL_SUMMATION_ORDER */

    // second leapfrog step as part of integration of equations of motion
    leapfrog_full_step(v, static_cast<vector_type>(f));

    // zero values for virtual particles to allow parallel reduction
    if (tag == VIRTUAL_PARTICLE) {
	v = 0;
	f = 0;
	en = 0;
	virial = 0;
    }

    // store particle associated with this thread
    g_v[GTID] = v;
    g_f[GTID] = static_cast<vector_type>(f);
    g_en[GTID] = en;
    g_virial[GTID] = virial;
}

/**
 * assign particles to cells
 */
template <typename T>
__global__ void assign_cells(float4 const* g_ir, float4* g_or, unsigned int* g_otag)
{
    __shared__ T s_ir[CELL_SIZE];
    __shared__ T s_or[CELL_SIZE];
    __shared__ unsigned int s_itag[CELL_SIZE];
    __shared__ unsigned int s_otag[CELL_SIZE];
    __shared__ unsigned int s_cell[CELL_SIZE];

    // number of particles in cell
    uint n = 0;
    // mark all particles in cell as virtual particles
    s_otag[threadIdx.x] = VIRTUAL_PARTICLE;
    __syncthreads();

    for (uint i = 0; i < npart; i += CELL_SIZE) {
	// load block of particles from global device memory
	T r = detach_particle_tag(g_ir[i + threadIdx.x], s_itag[threadIdx.x]);
	s_ir[threadIdx.x] = r;
	s_cell[threadIdx.x] = compute_cell(r);
	__syncthreads();

	if (threadIdx.x == 0) {
	    for (uint j = 0; j < CELL_SIZE && (i + j) < npart; j++) {
		if (s_cell[j] == BID) {
		    // store particle in cell
		    s_or[n] = s_ir[j];
		    // store particle tag
		    s_otag[n] = s_itag[j];
		    // increment particle count in cell
		    ++n;
		}
	    }
	}
	__syncthreads();
    }

    // store cell in global device memory
    unsigned int const tag = s_otag[threadIdx.x];
    g_or[BID * CELL_SIZE + threadIdx.x] = attach_particle_tag(s_or[threadIdx.x], tag);
    g_otag[BID * CELL_SIZE + threadIdx.x] = tag;
}

/**
 * examine neighbour cell for particles which moved into this block's cell
 */
template <typename T, typename U, typename I>
__device__ void examine_cell(I const& offset, float4 const* g_ir, U const* g_iR, U const* g_iv, T* s_or, T* s_oR, T* s_ov, unsigned int* s_otag, uint& npart)
{
    __shared__ T s_ir[CELL_SIZE];
    __shared__ T s_iR[CELL_SIZE];
    __shared__ T s_iv[CELL_SIZE];
    __shared__ unsigned int s_itag[CELL_SIZE];
    __shared__ uint s_cell[CELL_SIZE];

    // shared memory barrier
    __syncthreads();

    // compute cell index
    uint cell = compute_neighbour_cell(offset);
    // load particles in cell from global device memory
    T r = detach_particle_tag(g_ir[cell * CELL_SIZE + threadIdx.x], s_itag[threadIdx.x]);
    s_ir[threadIdx.x] = r;
    s_iR[threadIdx.x] = g_iR[cell * CELL_SIZE + threadIdx.x];
    s_iv[threadIdx.x] = g_iv[cell * CELL_SIZE + threadIdx.x];
    // compute new cell
    s_cell[threadIdx.x] = compute_cell(r);
    __syncthreads();

    if (threadIdx.x == 0) {
	for (uint j = 0; j < CELL_SIZE; j++) {
	    // skip virtual particles
	    if (s_itag[j] == VIRTUAL_PARTICLE) break;

	    // if particle belongs to this cell
	    if (s_cell[j] == BID && npart < CELL_SIZE) {
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
}

/**
 * update cells
 */
template <typename T, typename U>
__global__ void update_cells(float4 const* g_ir, U const* g_iR, U const* g_iv, float4* g_or, U* g_oR, U* g_ov, unsigned int* g_otag)
{
    enum { dimension = T::static_size };

    __shared__ T s_or[CELL_SIZE];
    __shared__ T s_oR[CELL_SIZE];
    __shared__ T s_ov[CELL_SIZE];
    __shared__ unsigned int s_otag[CELL_SIZE];
    // number of particles in cell
    uint n = 0;

    // mark all particles in cell as virtual particles
    s_otag[threadIdx.x] = VIRTUAL_PARTICLE;
    __syncthreads();

    if (dimension == 3) {
	examine_cell(make_int3( 0,  0,  0), g_ir, g_iR, g_iv, s_or, s_oR, s_ov, s_otag, n);
	// visit 26 neighbour cells
	examine_cell(make_int3(-1,  0,  0), g_ir, g_iR, g_iv, s_or, s_oR, s_ov, s_otag, n);
	examine_cell(make_int3(+1,  0,  0), g_ir, g_iR, g_iv, s_or, s_oR, s_ov, s_otag, n);
	examine_cell(make_int3( 0, -1,  0), g_ir, g_iR, g_iv, s_or, s_oR, s_ov, s_otag, n);
	examine_cell(make_int3( 0, +1,  0), g_ir, g_iR, g_iv, s_or, s_oR, s_ov, s_otag, n);
	examine_cell(make_int3(-1, -1,  0), g_ir, g_iR, g_iv, s_or, s_oR, s_ov, s_otag, n);
	examine_cell(make_int3(-1, +1,  0), g_ir, g_iR, g_iv, s_or, s_oR, s_ov, s_otag, n);
	examine_cell(make_int3(+1, -1,  0), g_ir, g_iR, g_iv, s_or, s_oR, s_ov, s_otag, n);
	examine_cell(make_int3(+1, +1,  0), g_ir, g_iR, g_iv, s_or, s_oR, s_ov, s_otag, n);
	examine_cell(make_int3( 0,  0, -1), g_ir, g_iR, g_iv, s_or, s_oR, s_ov, s_otag, n);
	examine_cell(make_int3(-1,  0, -1), g_ir, g_iR, g_iv, s_or, s_oR, s_ov, s_otag, n);
	examine_cell(make_int3(+1,  0, -1), g_ir, g_iR, g_iv, s_or, s_oR, s_ov, s_otag, n);
	examine_cell(make_int3( 0, -1, -1), g_ir, g_iR, g_iv, s_or, s_oR, s_ov, s_otag, n);
	examine_cell(make_int3( 0, +1, -1), g_ir, g_iR, g_iv, s_or, s_oR, s_ov, s_otag, n);
	examine_cell(make_int3(-1, -1, -1), g_ir, g_iR, g_iv, s_or, s_oR, s_ov, s_otag, n);
	examine_cell(make_int3(-1, +1, -1), g_ir, g_iR, g_iv, s_or, s_oR, s_ov, s_otag, n);
	examine_cell(make_int3(+1, -1, -1), g_ir, g_iR, g_iv, s_or, s_oR, s_ov, s_otag, n);
	examine_cell(make_int3(+1, +1, -1), g_ir, g_iR, g_iv, s_or, s_oR, s_ov, s_otag, n);
	examine_cell(make_int3( 0,  0, +1), g_ir, g_iR, g_iv, s_or, s_oR, s_ov, s_otag, n);
	examine_cell(make_int3(-1,  0, +1), g_ir, g_iR, g_iv, s_or, s_oR, s_ov, s_otag, n);
	examine_cell(make_int3(+1,  0, +1), g_ir, g_iR, g_iv, s_or, s_oR, s_ov, s_otag, n);
	examine_cell(make_int3( 0, -1, +1), g_ir, g_iR, g_iv, s_or, s_oR, s_ov, s_otag, n);
	examine_cell(make_int3( 0, +1, +1), g_ir, g_iR, g_iv, s_or, s_oR, s_ov, s_otag, n);
	examine_cell(make_int3(-1, -1, +1), g_ir, g_iR, g_iv, s_or, s_oR, s_ov, s_otag, n);
	examine_cell(make_int3(-1, +1, +1), g_ir, g_iR, g_iv, s_or, s_oR, s_ov, s_otag, n);
	examine_cell(make_int3(+1, -1, +1), g_ir, g_iR, g_iv, s_or, s_oR, s_ov, s_otag, n);
	examine_cell(make_int3(+1, +1, +1), g_ir, g_iR, g_iv, s_or, s_oR, s_ov, s_otag, n);
    }
    else {
	examine_cell(make_int2( 0,  0), g_ir, g_iR, g_iv, s_or, s_oR, s_ov, s_otag, n);
	// visit 8 neighbour cells
	examine_cell(make_int2(-1,  0), g_ir, g_iR, g_iv, s_or, s_oR, s_ov, s_otag, n);
	examine_cell(make_int2(+1,  0), g_ir, g_iR, g_iv, s_or, s_oR, s_ov, s_otag, n);
	examine_cell(make_int2( 0, -1), g_ir, g_iR, g_iv, s_or, s_oR, s_ov, s_otag, n);
	examine_cell(make_int2( 0, +1), g_ir, g_iR, g_iv, s_or, s_oR, s_ov, s_otag, n);
	examine_cell(make_int2(-1, -1), g_ir, g_iR, g_iv, s_or, s_oR, s_ov, s_otag, n);
	examine_cell(make_int2(-1, +1), g_ir, g_iR, g_iv, s_or, s_oR, s_ov, s_otag, n);
	examine_cell(make_int2(+1, -1), g_ir, g_iR, g_iv, s_or, s_oR, s_ov, s_otag, n);
	examine_cell(make_int2(+1, +1), g_ir, g_iR, g_iv, s_or, s_oR, s_ov, s_otag, n);
    }

    // store cell in global device memory
    unsigned int const tag = s_otag[threadIdx.x];
    g_or[BID * CELL_SIZE + threadIdx.x] = attach_particle_tag(s_or[threadIdx.x], tag);
    g_oR[BID * CELL_SIZE + threadIdx.x] = s_oR[threadIdx.x];
    g_ov[BID * CELL_SIZE + threadIdx.x] = s_ov[threadIdx.x];
    g_otag[BID * CELL_SIZE + threadIdx.x] = tag;
}

/**
 * first leapfrog step of integration of equations of motion
 */
template <int dimension, typename T>
__global__ void inteq(float4* g_r, T* g_R, T* g_v, T const* g_f)
{
    vector<float, dimension> r, dr, R, v, f;
    unsigned int tag;
    r = detach_particle_tag(g_r[GTID], tag);
    R = g_R[GTID];
    v = g_v[GTID];
    f = g_f[GTID];

    leapfrog_half_step(r, dr, R, v, f);

    g_r[GTID] = attach_particle_tag(r, tag);
    g_R[GTID] = R;
    g_v[GTID] = v;
}

}}} // namespace ljgpu::gpu::ljfluid

namespace ljgpu { namespace gpu
{

typedef ljfluid_base<ljfluid_impl_gpu_cell> _Base;
typedef ljfluid<ljfluid_impl_gpu_cell, 3> _3D;
typedef ljfluid<ljfluid_impl_gpu_cell, 2> _2D;

/**
 * device constant wrappers
 */
cuda::symbol<uint> _Base::ncell(cu::ljfluid::ncell);

/**
 * device function wrappers
 */
cuda::function<void (float4 const*, float4*, float4*, float*, float4*)>
    _3D::template variant<UNARY, C0POT>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 3>, UNARY, C0POT>);
cuda::function<void (float4 const*, float4*, float4*, float*, float4*)>
    _3D::template variant<UNARY, C2POT>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 3>, UNARY, C2POT>);
cuda::function<void (float4 const*, float4*, float4*, float*, float4*)>
    _3D::template variant<BINARY, C0POT>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 3>, BINARY, C0POT>);
cuda::function<void (float4 const*, float4*, float4*, float*, float4*)>
    _3D::template variant<BINARY, C2POT>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 3>, BINARY, C2POT>);

cuda::function<void (float4 const*, float4*, unsigned int*)>
    _3D::assign_cells(cu::ljfluid::assign_cells<cu::vector<float, 3> >);
cuda::function<void (float4 const*, float4 const*, float4 const*, float4*, float4*, float4*, unsigned int*)>
    _3D::update_cells(cu::ljfluid::update_cells<cu::vector<float, 3> >);
cuda::function<void (float4*, float4*, float4*, float4 const*)>
    _3D::inteq(cu::ljfluid::inteq<3>);

cuda::function<void (float4 const*, float2*, float2*, float*, float2*)>
    _2D::template variant<UNARY, C0POT>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 2>, UNARY, C0POT>);
cuda::function<void (float4 const*, float2*, float2*, float*, float2*)>
    _2D::template variant<UNARY, C2POT>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 2>, UNARY, C2POT>);
cuda::function<void (float4 const*, float2*, float2*, float*, float2*)>
    _2D::template variant<BINARY, C0POT>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 2>, BINARY, C0POT>);
cuda::function<void (float4 const*, float2*, float2*, float*, float2*)>
    _2D::template variant<BINARY, C2POT>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 2>, BINARY, C2POT>);

cuda::function<void (float4 const*, float4*, unsigned int*)>
    _2D::assign_cells(cu::ljfluid::assign_cells<cu::vector<float, 2> >);
cuda::function<void (float4 const*, float2 const*, float2 const*, float4*, float2*, float2*, unsigned int*)>
    _2D::update_cells(cu::ljfluid::update_cells<cu::vector<float, 2> >);
cuda::function<void (float4*, float2*, float2*, float2 const*)>
    _2D::inteq(cu::ljfluid::inteq<2>);

}} // namespace ljgpu::gpu
