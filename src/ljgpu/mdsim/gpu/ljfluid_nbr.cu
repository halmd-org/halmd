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

/** number of cells per dimension */
__constant__ uint ncell;
/** neighbour list length */
__constant__ uint nbl_size;
/** neighbour list stride */
__constant__ uint nbl_stride;
/** squared potential cutoff distance with neighbour list skin */
__constant__ float rr_nbl;
/** neighbour lists in global device memory */
__constant__ unsigned int* g_nbl;

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
	  typename T>
__global__ void mdstep(float4 const* g_r, T* g_v, T* g_f, float* g_en, T* g_virial)
{
    enum { dimension = vector_type::static_size };

    // load particle associated with this thread
    unsigned int tag;
    vector_type r = detach_particle_tag(g_r[GTID], tag);
    // particle type in binary mixture
    int const a = (tag >= mpart[0]);

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

    for (uint i = 0; i < nbl_size; ++i) {
	// coalesced read from neighbour list
	unsigned int const n = g_nbl[i * nbl_stride + GTID];
	// skip placeholder particles
	if (n == VIRTUAL_PARTICLE) break;

	unsigned int tag_;
	vector_type r_ = detach_particle_tag(tex1Dfetch(tex<dimension>::r, n), tag_);
	// particle type in binary mixture
	int const b = (tag_ >= mpart[0]);

	// accumulate force between particles
	compute_force<mixture, potential>(r, r_, f, en, virial, a + b);
    }

#ifdef USE_VERLET_DSFUN
    vector<dsfloat, dimension> v(g_v[GTID], g_v[GTID + GTDIM]);
#else
    vector_type v = g_v[GTID];
#endif
    // second leapfrog step as part of integration of equations of motion
    leapfrog_full_step(v, static_cast<vector_type>(f));

    // store particle associated with this thread
    g_v[GTID] = static_cast<vector_type>(v);
#ifdef USE_VERLET_DSFUN
    g_v[GTID + GTDIM] = dsfloat2lo(v);
#endif
    g_f[GTID] = static_cast<vector_type>(f);
    g_en[GTID] = en;
    g_virial[GTID] = virial;
}

/**
 * compute neighbour cell
 */
__device__ uint compute_neighbour_cell(int3 const &offset)
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
 * update neighbour list with particles of given cell
 */
template <bool same_cell, typename T, typename I>
__device__ void update_cell_neighbours(I const& offset, unsigned int const* g_cell, T const& r, unsigned int const& n, uint& count)
{
    extern __shared__ unsigned int s_n[];
    T* const s_r = reinterpret_cast<T*>(&s_n[blockDim.x]);
    enum { dimension = T::static_size };

    // shared memory barrier
    __syncthreads();

    // compute cell index
    uint const cell = compute_neighbour_cell(offset);
    // load particles in cell
    unsigned int const n_ = g_cell[cell * blockDim.x + threadIdx.x];
    s_n[threadIdx.x] = n_;
    s_r[threadIdx.x] = tex1Dfetch(tex<dimension>::r, n_);
    __syncthreads();

    if (n == VIRTUAL_PARTICLE) return;

    for (uint i = 0; i < blockDim.x; ++i) {
	// particle number of cell placeholder
	unsigned int const m = s_n[i];
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
__global__ void update_neighbours(unsigned int* g_ret, unsigned int const* g_cell)
{
    // load particle from cell placeholder
    unsigned int const n = g_cell[GTID];
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

    // return failure if any neighbour list is fully occupied
    if (count == nbl_size) {
	*g_ret = EXIT_FAILURE;
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
__global__ void find_cell_offset(uint* g_cell, unsigned int* g_cell_offset)
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
__global__ void assign_cells(unsigned int* g_ret, uint const* g_cell, unsigned int const* g_cell_offset, unsigned int const* g_itag, unsigned int* g_otag)
{
    __shared__ unsigned int s_offset[1];

    if (threadIdx.x == 0) {
	s_offset[0] = g_cell_offset[BID];
    }
    __syncthreads();
    // global offset of first particle in this block's cell
    const unsigned int offset = s_offset[0];
    // global offset of this thread's particle
    const unsigned int n = offset + threadIdx.x;
    // mark as virtual particle
    unsigned int tag = VIRTUAL_PARTICLE;
    // mark as real particle if appropriate
    if (offset != VIRTUAL_PARTICLE && n < npart && g_cell[n] == BID) {
	tag = g_itag[n];
    }
    // return failure if any cell list is fully occupied
    if (tag != VIRTUAL_PARTICLE && (threadIdx.x + 1) == blockDim.x) {
	*g_ret = EXIT_FAILURE;
    }
    // store particle in this block's cell
    g_otag[BID * blockDim.x + threadIdx.x] = tag;
}

/**
 * generate ascending index sequence
 */
__global__ void gen_index(unsigned int* g_index)
{
    g_index[GTID] = (GTID < npart) ? GTID : 0;
}

/**
 * order particles after given permutation
 */
template <int dimension, typename T>
__global__ void order_particles(unsigned int const* g_index, float4* g_or, T* g_oR, T* g_ov, unsigned int* g_otag)
{
    // permutation index
    uint const j = g_index[GTID];
    // permute particle phase space coordinates
    unsigned int tag;
    vector<float, dimension> r = detach_particle_tag(tex1Dfetch(tex<dimension>::r, j), tag);
    g_or[GTID] = attach_particle_tag(r, tag);
    g_oR[GTID] = tex1Dfetch(tex<dimension>::R, j);
    g_ov[GTID] = tex1Dfetch(tex<dimension>::v, j);
#ifdef USE_VERLET_DSFUN
    g_or[GTID + GTDIM] = tex1Dfetch(tex<dimension>::r, j + GTDIM);
    g_ov[GTID + GTDIM] = tex1Dfetch(tex<dimension>::v, j + GTDIM);
#endif
    g_otag[GTID] = tag;
}

/**
 * assign velocities given particle order in memory
 */
template <int dimension, typename T>
__global__ void order_velocities(uint const* g_tag, T* g_v)
{
    g_v[GTID] = tex1Dfetch(tex<dimension>::v, g_tag[GTID]);
}

/**
 * sample trajectories
 */
template <int dimension, typename T>
__global__ void sample(unsigned int const* g_perm, T* g_or, T* g_ov)
{
    // permutation index
    uint const j = g_perm[GTID];
    // permute particle phase space coordinates
    vector<float, dimension> const r = tex1Dfetch(tex<dimension>::r, j);
    vector<float, dimension> const R = tex1Dfetch(tex<dimension>::R, j);
    g_or[GTID] = r + box * R;
    g_ov[GTID] = tex1Dfetch(tex<dimension>::v, j);
}

/**
 * first leapfrog step of integration of equations of motion
 */
template <int dimension, typename T>
__global__ void inteq(float4* g_r, T* g_dr, T* g_R, T* g_v, T const* g_f)
{
    unsigned int tag;
#ifdef USE_VERLET_DSFUN
    vector<dsfloat, dimension> r(detach_particle_tag(g_r[GTID], tag), g_r[GTID + GTDIM]);
    vector<dsfloat, dimension> v(g_v[GTID], g_v[GTID + GTDIM]);
    vector<dsfloat, dimension> dr;
#else
    vector<float, dimension> r = detach_particle_tag(g_r[GTID], tag);
    vector<float, dimension> v = g_v[GTID];
    vector<float, dimension> dr;
#endif
    vector<float, dimension> R = g_R[GTID];
    vector<float, dimension> f = g_f[GTID];

    leapfrog_half_step(r, dr, R, v, f);

#ifdef USE_VERLET_DSFUN
    g_r[GTID] = attach_particle_tag(dsfloat2hi(r), tag);
    g_r[GTID + GTDIM] = dsfloat2lo(r);
    // particle displacement for neighbour list update constraint
    g_dr[GTID] = dsfloat2hi(dr) + g_dr[GTID];
    g_v[GTID] = dsfloat2hi(v);
    g_v[GTID + GTDIM] = dsfloat2lo(v);
#else
    g_r[GTID] = attach_particle_tag(r, tag);
    g_dr[GTID] = dr + g_dr[GTID];
    g_v[GTID] = v;
#endif
    g_R[GTID] = R;
}

}}} // namespace ljgpu::gpu::ljfluid

namespace ljgpu { namespace gpu
{

typedef ljfluid_base<ljfluid_impl_gpu_neighbour> _Base;
typedef ljfluid<ljfluid_impl_gpu_neighbour, 3> _3D;
typedef ljfluid<ljfluid_impl_gpu_neighbour, 2> _2D;

/**
 * device constant wrappers
 */
cuda::symbol<uint> _Base::ncell(cu::ljfluid::ncell);
cuda::symbol<uint> _Base::nbl_size(cu::ljfluid::nbl_size);
cuda::symbol<uint> _Base::nbl_stride(cu::ljfluid::nbl_stride);
cuda::symbol<float> _Base::rr_nbl(cu::ljfluid::rr_nbl);
cuda::symbol<unsigned int*> _Base::g_nbl(cu::ljfluid::g_nbl);

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
cuda::function<void (unsigned int*, uint const*, unsigned int const*, unsigned int const*, unsigned int*)>
    _Base::assign_cells(cu::ljfluid::assign_cells);
cuda::function<void (uint*, unsigned int*)>
    _Base::find_cell_offset(cu::ljfluid::find_cell_offset);
cuda::function<void (unsigned int*)>
    _Base::gen_index(cu::ljfluid::gen_index);

cuda::function<void (float4 const*, float4*, float4*, float*, float4*)>
    _3D::template variant<UNARY, C0POT>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 3>, UNARY, C0POT>);
cuda::function<void (float4 const*, float4*, float4*, float*, float4*)>
    _3D::template variant<UNARY, C2POT>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 3>, UNARY, C2POT>);
cuda::function<void (float4 const*, float4*, float4*, float*, float4*)>
    _3D::template variant<BINARY, C0POT>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 3>, BINARY, C0POT>);
cuda::function<void (float4 const*, float4*, float4*, float*, float4*)>
    _3D::template variant<BINARY, C2POT>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 3>, BINARY, C2POT>);

cuda::function<void (unsigned int*, unsigned int const*)>
    _3D::update_neighbours(cu::ljfluid::update_neighbours<3>);
cuda::function<void (float4 const*, uint*)>
    _3D::compute_cell(cu::ljfluid::compute_cell<3>);
cuda::function<void (unsigned int const*, float4*, float4*, float4*, unsigned int*)>
    _3D::order_particles(cu::ljfluid::order_particles<3>);
cuda::function<void (uint const*, float4*)>
    _3D::order_velocities(cu::ljfluid::order_velocities<3>);
cuda::function<void (unsigned int const*, float4*, float4*)>
    _3D::sample(cu::ljfluid::sample<3>);
cuda::function<void (float4*, float4*, float4*, float4*, float4 const*)>
    _3D::inteq(cu::ljfluid::inteq<3>);

cuda::function<void (float4 const*, float2*, float2*, float*, float2*)>
    _2D::template variant<UNARY, C0POT>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 2>, UNARY, C0POT>);
cuda::function<void (float4 const*, float2*, float2*, float*, float2*)>
    _2D::template variant<UNARY, C2POT>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 2>, UNARY, C2POT>);
cuda::function<void (float4 const*, float2*, float2*, float*, float2*)>
    _2D::template variant<BINARY, C0POT>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 2>, BINARY, C0POT>);
cuda::function<void (float4 const*, float2*, float2*, float*, float2*)>
    _2D::template variant<BINARY, C2POT>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 2>, BINARY, C2POT>);

cuda::function<void (unsigned int*, unsigned int const*)>
    _2D::update_neighbours(cu::ljfluid::update_neighbours<2>);
cuda::function<void (float4 const*, uint*)>
    _2D::compute_cell(cu::ljfluid::compute_cell<2>);
cuda::function<void (unsigned int const*, float4*, float2*, float2*, unsigned int*)>
    _2D::order_particles(cu::ljfluid::order_particles<2>);
cuda::function<void (uint const*, float2*)>
    _2D::order_velocities(cu::ljfluid::order_velocities<2>);
cuda::function<void (unsigned int const*, float2*, float2*)>
    _2D::sample(cu::ljfluid::sample<2>);
cuda::function<void (float4*, float2*, float2*, float2*, float2 const*)>
    _2D::inteq(cu::ljfluid::inteq<2>);

}} // namespace ljgpu::gpu
