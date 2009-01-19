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
#include <ljgpu/mdsim/gpu/ljfluid_nbr.hpp>
using namespace ljgpu::gpu;

namespace ljgpu { namespace cu { namespace ljfluid
{

/** fixed number of placeholders per cell */
enum { CELL_SIZE = ljfluid_base<ljfluid_impl_gpu_neighbour>::CELL_SIZE };

/** number of cells per dimension */
__constant__ uint ncell;
/** neighbour list length */
__constant__ uint nbl_size;
/** neighbour list stride */
__constant__ uint nbl_stride;
/** squared potential cutoff distance with neighbour list skin */
__constant__ float rr_nbl;
/** neighbour lists in global device memory */
__constant__ int* g_nbl;

/** n-dimensional particle texture references */
template <int dimension>
struct tex;

template <>
struct tex<3>
{
    /** periodic particle positions */
    static texture<float4, 1, cudaReadModeElementType> r;
    /** extended particle positions */
    static texture<float4, 1, cudaReadModeElementType> R;
    /** particle velocities */
    static texture<float4, 1, cudaReadModeElementType> v;
};

template <>
struct tex<2>
{
    /** periodic particle positions */
    static texture<float4, 1, cudaReadModeElementType> r;
    /** extended particle positions */
    static texture<float2, 1, cudaReadModeElementType> R;
    /** particle velocities */
    static texture<float2, 1, cudaReadModeElementType> v;
};

// instantiate texture references
texture<float4, 1, cudaReadModeElementType> tex<3>::r;
texture<float4, 1, cudaReadModeElementType> tex<2>::r;
texture<float4, 1, cudaReadModeElementType> tex<3>::R;
texture<float2, 1, cudaReadModeElementType> tex<2>::R;
texture<float4, 1, cudaReadModeElementType> tex<3>::v;
texture<float2, 1, cudaReadModeElementType> tex<2>::v;

/**
 * n-dimensional MD simulation step
 */
template <typename vector_type,
          mixture_type mixture,
	  potential_type potential,
	  ensemble_type ensemble,
	  typename T>
__global__ void mdstep(float4 const* g_r, T* g_v, T* g_f, float* g_en, float* g_virial)
{
    enum { dimension = vector_type::static_size };

    // load particle associated with this thread
    vector_type r, v;
    int tag;
    (r, tag) = g_r[GTID];
    v = g_v[GTID];
    // extract particle type from tag
    uint const a = (tag >> 31);

    // potential energy contribution
    float en = 0;
    // virial equation sum contribution
    float virial = 0;
    // force sum
    vector<dfloat, dimension> f = 0;

    for (uint i = 0; i < nbl_size; ++i) {
	// coalesced read from neighbour list
	int const n = g_nbl[i * nbl_stride + GTID];
	// skip placeholder particles
	if (n == VIRTUAL_PARTICLE) break;

	vector_type r_;
	int tag_;
	(r_, tag_) = tex1Dfetch(tex<dimension>::r, n);
	// particle type
	uint const b = (tag_ >> 31);

	// accumulate force between particles
	compute_force<mixture, potential>(r, r_, f, en, virial, a + b);
    }

    // second leapfrog step as part of integration of equations of motion
    leapfrog_full_step(v, static_cast<vector_type>(f));

    // random collisions with heat bath
    if (ensemble == NVT) {
	anderson_thermostat(v);
    }

    // store particle associated with this thread
    g_v[GTID] = v;
    g_f[GTID] = static_cast<vector_type>(f);
    g_en[GTID] = en;
    g_virial[GTID] = virial;
}

/**
 * initialise particle tags
 */
template <typename vector_type>
__global__ void init_tags(float4* g_r, int* g_tag)
{
    vector_type const r = g_r[GTID];
    int tag = VIRTUAL_PARTICLE;
    if (GTID < npart) {
	tag = GTID;
    }
    g_r[GTID] = (r, tag);
    g_tag[GTID] = tag;
}

/**
 * compute neighbour cell
 */
__device__ uint compute_neighbour_cell(int3 const &offset)
{
    // cell belonging to this execution block
    int x = blockIdx.x % ncell;
    int y = (blockIdx.x / ncell) % ncell;
    int z = blockIdx.x / ncell / ncell;
    // neighbour cell of this cell
    x = (x + ncell + offset.x) % ncell;
    y = (y + ncell + offset.y) % ncell;
    z = (z + ncell + offset.z) % ncell;

    return (z * ncell + y) * ncell + x;
}

__device__ uint compute_neighbour_cell(int2 const& offset)
{
    // cell belonging to this execution block
    int x = blockIdx.x % ncell;
    int y = blockIdx.x / ncell;
    // neighbour cell of this cell
    x = (x + ncell + offset.x) % ncell;
    y = (y + ncell + offset.y) % ncell;

    return y * ncell + x;
}

/**
 * update neighbour list with particles of given cell
 */
template <bool same_cell, typename T, typename I>
__device__ void update_cell_neighbours(I const& offset, int const* g_cell, T const& r, int const& n, uint& count)
{
    __shared__ int s_n[CELL_SIZE];
    __shared__ T s_r[CELL_SIZE];
    enum { dimension = T::static_size };

    // shared memory barrier
    __syncthreads();

    // compute cell index
    uint const cell = compute_neighbour_cell(offset);
    // load particles in cell
    int const n_ = g_cell[cell * CELL_SIZE + threadIdx.x];
    s_n[threadIdx.x] = n_;
    s_r[threadIdx.x] = tex1Dfetch(tex<dimension>::r, n_);
    __syncthreads();

    if (n == VIRTUAL_PARTICLE) return;

    for (uint i = 0; i < CELL_SIZE; ++i) {
	// particle number of cell placeholder
	int const m = s_n[i];
	// skip placeholder particles
	if (m == VIRTUAL_PARTICLE) break;
	// skip same particle
	if (same_cell && i == threadIdx.x) continue;

	// particle distance vector
	T dr = r - s_r[i];
	// enforce periodic boundary conditions
	dr -= rintf(__fdividef(dr, box)) * box;
	// squared particle distance
	float rr = dr * dr;

	// enforce cutoff length with neighbour list skin
	if (rr <= rr_nbl && count < nbl_size) {
	    // scattered write to neighbour list
	    g_nbl[count * nbl_stride + n] = m;
	    // increment neighbour list particle count
	    count++;
	}
    }
}

/**
 * update neighbour lists
 */
template <uint dimension>
__global__ void update_neighbours(int const* g_cell)
{
    // load particle from cell placeholder
    int const n = g_cell[GTID];
    vector<float, dimension> const r = tex1Dfetch(tex<dimension>::r, n);
    // number of particles in neighbour list
    uint count = 0;

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
	// visit this cell
	update_cell_neighbours<true>(make_int3( 0,  0,  0), g_cell, r, n, count);
	// visit 26 neighbour cells, grouped into 13 pairs of mutually opposite cells
	update_cell_neighbours<false>(make_int3(-1, -1, -1), g_cell, r, n, count);
	update_cell_neighbours<false>(make_int3(+1, +1, +1), g_cell, r, n, count);
	update_cell_neighbours<false>(make_int3(-1, -1, +1), g_cell, r, n, count);
	update_cell_neighbours<false>(make_int3(+1, +1, -1), g_cell, r, n, count);
	update_cell_neighbours<false>(make_int3(-1, +1, +1), g_cell, r, n, count);
	update_cell_neighbours<false>(make_int3(+1, -1, -1), g_cell, r, n, count);
	update_cell_neighbours<false>(make_int3(+1, -1, +1), g_cell, r, n, count);
	update_cell_neighbours<false>(make_int3(-1, +1, -1), g_cell, r, n, count);
	update_cell_neighbours<false>(make_int3(-1, -1,  0), g_cell, r, n, count);
	update_cell_neighbours<false>(make_int3(+1, +1,  0), g_cell, r, n, count);
	update_cell_neighbours<false>(make_int3(-1, +1,  0), g_cell, r, n, count);
	update_cell_neighbours<false>(make_int3(+1, -1,  0), g_cell, r, n, count);
	update_cell_neighbours<false>(make_int3(-1,  0, -1), g_cell, r, n, count);
	update_cell_neighbours<false>(make_int3(+1,  0, +1), g_cell, r, n, count);
	update_cell_neighbours<false>(make_int3(-1,  0, +1), g_cell, r, n, count);
	update_cell_neighbours<false>(make_int3(+1,  0, -1), g_cell, r, n, count);
	update_cell_neighbours<false>(make_int3( 0, -1, -1), g_cell, r, n, count);
	update_cell_neighbours<false>(make_int3( 0, +1, +1), g_cell, r, n, count);
	update_cell_neighbours<false>(make_int3( 0, -1, +1), g_cell, r, n, count);
	update_cell_neighbours<false>(make_int3( 0, +1, -1), g_cell, r, n, count);
	update_cell_neighbours<false>(make_int3(-1,  0,  0), g_cell, r, n, count);
	update_cell_neighbours<false>(make_int3(+1,  0,  0), g_cell, r, n, count);
	update_cell_neighbours<false>(make_int3( 0, -1,  0), g_cell, r, n, count);
	update_cell_neighbours<false>(make_int3( 0, +1,  0), g_cell, r, n, count);
	update_cell_neighbours<false>(make_int3( 0,  0, -1), g_cell, r, n, count);
	update_cell_neighbours<false>(make_int3( 0,  0, +1), g_cell, r, n, count);
    }
    else {
	// visit this cell
	update_cell_neighbours<true>(make_int2( 0,  0), g_cell, r, n, count);
	// visit 8 neighbour cells, grouped into 4 pairs of mutually opposite cells
	update_cell_neighbours<false>(make_int2(-1, -1), g_cell, r, n, count);
	update_cell_neighbours<false>(make_int2(+1, +1), g_cell, r, n, count);
	update_cell_neighbours<false>(make_int2(-1, +1), g_cell, r, n, count);
	update_cell_neighbours<false>(make_int2(+1, -1), g_cell, r, n, count);
	update_cell_neighbours<false>(make_int2(-1,  0), g_cell, r, n, count);
	update_cell_neighbours<false>(make_int2(+1,  0), g_cell, r, n, count);
	update_cell_neighbours<false>(make_int2( 0, -1), g_cell, r, n, count);
	update_cell_neighbours<false>(make_int2( 0, +1), g_cell, r, n, count);
    }
}

/**
 * compute cell indices for given particle positions
 */
__device__ uint compute_cell_index(vector<float, 3> r)
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
    return uint(r.x) + ncell * (uint(r.y) + ncell * uint(r.z));
}

__device__ uint compute_cell_index(vector<float, 2> r)
{
    r = (__saturatef(r / box) * (1.f - FLT_EPSILON)) * ncell;
    return uint(r.x) + ncell * uint(r.y);
}

template <uint dimension>
__global__ void compute_cell(float4 const* g_part, uint* g_cell)
{
    vector<float, dimension> const r = g_part[GTID];
    g_cell[GTID] = compute_cell_index(r);
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
__global__ void assign_cells(uint const* g_cell, int const* g_cell_offset, int const* g_itag, int* g_otag)
{
    __shared__ int s_offset[1];

    if (threadIdx.x == 0) {
	s_offset[0] = g_cell_offset[blockIdx.x];
    }
    __syncthreads();
    // global offset of first particle in this block's cell
    const int offset = s_offset[0];
    // global offset of this thread's particle
    const int n = offset + threadIdx.x;
    // mark as virtual particle
    int tag = -1;
    // mark as real particle if appropriate
    if (offset >= 0 && n < npart && g_cell[n] == blockIdx.x) {
	tag = g_itag[n];
    }
    // store particle in this block's cell
    g_otag[blockIdx.x * CELL_SIZE + threadIdx.x] = tag;
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
template <int dimension, typename T>
__global__ void order_particles(const int* g_idx, float4* g_or, T* g_oR, T* g_ov, int* g_tag)
{
    // permutation index
    uint const j = g_idx[GTID];
    // permute particle phase space coordinates
    vector<float, dimension> r;
    int tag;
    (r, tag) = tex1Dfetch(tex<dimension>::r, j);
    g_or[GTID] = (r, tag);
    g_oR[GTID] = tex1Dfetch(tex<dimension>::R, j);
    g_ov[GTID] = tex1Dfetch(tex<dimension>::v, j);
    g_tag[GTID] = tag;
}

}}} // namespace ljgpu::gpu::ljfluid

namespace ljgpu { namespace gpu
{

typedef ljfluid_base<ljfluid_impl_gpu_neighbour> _Base;
typedef ljfluid<ljfluid_impl_gpu_neighbour<3> > _3D;
typedef ljfluid<ljfluid_impl_gpu_neighbour<2> > _2D;

/**
 * device constant wrappers
 */
cuda::symbol<uint> _Base::ncell(cu::ljfluid::ncell);
cuda::symbol<uint> _Base::nbl_size(cu::ljfluid::nbl_size);
cuda::symbol<uint> _Base::nbl_stride(cu::ljfluid::nbl_stride);
cuda::symbol<float> _Base::rr_nbl(cu::ljfluid::rr_nbl);
cuda::symbol<int*> _Base::g_nbl(cu::ljfluid::g_nbl);

/**
 * device texture wrappers
 */
cuda::texture<float4> _3D::r(cu::ljfluid::tex<3>::r);
cuda::texture<float4> _2D::r(cu::ljfluid::tex<2>::r);
cuda::texture<float4> _3D::R(cu::ljfluid::tex<3>::R);
cuda::texture<float2> _2D::R(cu::ljfluid::tex<2>::R);
cuda::texture<float4> _3D::v(cu::ljfluid::tex<3>::v);
cuda::texture<float2> _2D::v(cu::ljfluid::tex<2>::v);

/**
 * device function wrappers
 */
cuda::function<void (uint const*, int const*, int const*, int*)>
    _Base::assign_cells(cu::ljfluid::assign_cells);
cuda::function<void (uint*, int*)>
    _Base::find_cell_offset(cu::ljfluid::find_cell_offset);
cuda::function<void (int*)>
    _Base::gen_index(cu::ljfluid::gen_index);

cuda::function<void (float4 const*, float4*, float4*, float*, float*)>
    _3D::template variant<UNARY, C0POT, NVE>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 3>, UNARY, C0POT, NVE>);
cuda::function<void (float4 const*, float4*, float4*, float*, float*)>
    _3D::template variant<UNARY, C0POT, NVT>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 3>, UNARY, C0POT, NVT>);
cuda::function<void (float4 const*, float4*, float4*, float*, float*)>
    _3D::template variant<UNARY, C2POT, NVE>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 3>, UNARY, C2POT, NVE>);
cuda::function<void (float4 const*, float4*, float4*, float*, float*)>
    _3D::template variant<UNARY, C2POT, NVT>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 3>, UNARY, C2POT, NVT>);

cuda::function<void (float4 const*, float4*, float4*, float*, float*)>
    _3D::template variant<BINARY, C0POT, NVE>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 3>, BINARY, C0POT, NVE>);
cuda::function<void (float4 const*, float4*, float4*, float*, float*)>
    _3D::template variant<BINARY, C0POT, NVT>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 3>, BINARY, C0POT, NVT>);
cuda::function<void (float4 const*, float4*, float4*, float*, float*)>
    _3D::template variant<BINARY, C2POT, NVE>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 3>, BINARY, C2POT, NVE>);
cuda::function<void (float4 const*, float4*, float4*, float*, float*)>
    _3D::template variant<BINARY, C2POT, NVT>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 3>, BINARY, C2POT, NVT>);

cuda::function<void (float4*, int*)>
    _3D::init_tags(cu::ljfluid::init_tags<cu::vector<float, 3> >);
cuda::function<void (int const*)>
    _3D::update_neighbours(cu::ljfluid::update_neighbours<3>);
cuda::function<void (float4 const*, uint*)>
    _3D::compute_cell(cu::ljfluid::compute_cell<3>);
cuda::function<void (const int*, float4*, float4*, float4*, int*)>
    _3D::order_particles(cu::ljfluid::order_particles<3>);

cuda::function<void (float4 const*, float2*, float2*, float*, float*)>
    _2D::template variant<UNARY, C0POT, NVE>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 2>, UNARY, C0POT, NVE>);
cuda::function<void (float4 const*, float2*, float2*, float*, float*)>
    _2D::template variant<UNARY, C0POT, NVT>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 2>, UNARY, C0POT, NVT>);
cuda::function<void (float4 const*, float2*, float2*, float*, float*)>
    _2D::template variant<UNARY, C2POT, NVE>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 2>, UNARY, C2POT, NVE>);
cuda::function<void (float4 const*, float2*, float2*, float*, float*)>
    _2D::template variant<UNARY, C2POT, NVT>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 2>, UNARY, C2POT, NVT>);

cuda::function<void (float4 const*, float2*, float2*, float*, float*)>
    _2D::template variant<BINARY, C0POT, NVE>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 2>, BINARY, C0POT, NVE>);
cuda::function<void (float4 const*, float2*, float2*, float*, float*)>
    _2D::template variant<BINARY, C0POT, NVT>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 2>, BINARY, C0POT, NVT>);
cuda::function<void (float4 const*, float2*, float2*, float*, float*)>
    _2D::template variant<BINARY, C2POT, NVE>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 2>, BINARY, C2POT, NVE>);
cuda::function<void (float4 const*, float2*, float2*, float*, float*)>
    _2D::template variant<BINARY, C2POT, NVT>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 2>, BINARY, C2POT, NVT>);

cuda::function<void (float4*, int*)>
    _2D::init_tags(cu::ljfluid::init_tags<cu::vector<float, 2> >);
cuda::function<void (int const*)>
    _2D::update_neighbours(cu::ljfluid::update_neighbours<2>);
cuda::function<void (float4 const*, uint*)>
    _2D::compute_cell(cu::ljfluid::compute_cell<2>);
cuda::function<void (const int*, float4*, float2*, float2*, int*)>
    _2D::order_particles(cu::ljfluid::order_particles<2>);

}} // namespace ljgpu::gpu
