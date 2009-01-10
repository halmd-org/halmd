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
#include <ljgpu/ljfluid/gpu/base.cuh>
#include <ljgpu/ljfluid/gpu/ljfluid_nbr.hpp>

namespace ljgpu { namespace gpu
{
enum {
    /** fixed number of placeholders per cell */
    CELL_SIZE = ljfluid_base<ljfluid_impl_gpu_neighbour>::CELL_SIZE,
    /** virtual particle tag */
    VIRTUAL_PARTICLE = ljfluid_base<ljfluid_impl_gpu_neighbour>::VIRTUAL_PARTICLE,
};


/** number of cells per dimension */
static __constant__ uint ncell;
/** neighbour list length */
static __constant__ uint nbl_size;
/** neighbour list stride */
static __constant__ uint nbl_stride;
/** potential cutoff radius with cell skin */
static __constant__ float r_cell;
/** squared potential radius distance with cell skin */
static __constant__ float rr_cell;

/** n-dimensional particle texture references */
template <typename T = void>
struct tex;

template <>
struct tex<void>
{
    /** texture reference to particle tags */
    static texture<int, 1, cudaReadModeElementType> tag;
};

template <>
struct tex<float4>
{
    /** periodic particle positions */
    static texture<float4, 1, cudaReadModeElementType> r;
    /** extended particle positions */
    static texture<float4, 1, cudaReadModeElementType> R;
    /** particle velocities */
    static texture<float4, 1, cudaReadModeElementType> v;
};

template <>
struct tex<float2>
{
    /** periodic particle positions */
    static texture<float2, 1, cudaReadModeElementType> r;
    /** extended particle positions */
    static texture<float2, 1, cudaReadModeElementType> R;
    /** particle velocities */
    static texture<float2, 1, cudaReadModeElementType> v;
};

// instantiate texture references
texture<int, 1, cudaReadModeElementType> tex<>::tag;
texture<float4, 1, cudaReadModeElementType> tex<float4>::r;
texture<float4, 1, cudaReadModeElementType> tex<float4>::R;
texture<float4, 1, cudaReadModeElementType> tex<float4>::v;
texture<float2, 1, cudaReadModeElementType> tex<float2>::r;
texture<float2, 1, cudaReadModeElementType> tex<float2>::R;
texture<float2, 1, cudaReadModeElementType> tex<float2>::v;

/**
 * n-dimensional MD simulation step
 */
template <typename T, typename TT, typename U>
__global__ void mdstep(U const* g_r, U* g_v, U* g_f, int const* g_nbl, float* g_en, float* g_virial)
{
    // load particle associated with this thread
    const T r = unpack(g_r[GTID]);
    T v = unpack(g_v[GTID]);

    // potential energy contribution
    float en = 0;
    // virial equation sum contribution
    float virial = 0;
    // force sum
    TT f = 0;

    for (uint i = 0; i < nbl_size; ++i) {
	// coalesced read from neighbour list
	const int n = g_nbl[i * nbl_stride + GTID];
	// skip placeholder particles
	if (n == VIRTUAL_PARTICLE)
	    break;
	// accumulate force between particles
	compute_force(r, unpack(tex1Dfetch(tex<U>::r, n)), f, en, virial);
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
 * initialise particle tags
 */
static __global__ void init_tags(int* g_tag)
{
    int tag = VIRTUAL_PARTICLE;
    if (GTID < npart) {
	tag = GTID;
    }
    g_tag[GTID] = tag;
}

/**
 * compute neighbour cell
 */
__device__ uint compute_neighbour_cell(int3 const &offset)
{
    int3 i = make_int3(blockIdx.x % ncell,
		       (blockIdx.x / ncell) % ncell,
		       blockIdx.x / ncell / ncell);
    i = make_int3((i.x + ncell + offset.x) % ncell,
		  (i.y + ncell + offset.y) % ncell,
		  (i.z + ncell + offset.z) % ncell);
    return (i.z * ncell + i.y) * ncell + i.x;
}

__device__ uint compute_neighbour_cell(int2 const& offset)
{
    int2 i = make_int2(blockIdx.x % ncell,
		       blockIdx.x / ncell);
    i = make_int2((i.x + ncell + offset.x) % ncell,
		  (i.y + ncell + offset.y) % ncell);
    return i.y * ncell + i.x;
}

/**
 * update neighbour list with particles of given cell
 */
template <uint cell_size, typename U, bool same_cell, typename T, typename I>
__device__ void update_cell_neighbours(I const& offset, int const* g_cell, int* g_nbl, T const& r, int const& n, uint& count)
{
    __shared__ int s_n[cell_size];
    __shared__ T s_r[cell_size];

    // compute cell index
    const uint cell = compute_neighbour_cell(offset);

    // load particles in cell
    s_n[threadIdx.x] = g_cell[cell * cell_size + threadIdx.x];
    s_r[threadIdx.x] = unpack(tex1Dfetch(tex<U>::r, s_n[threadIdx.x]));
    __syncthreads();

    if (n != VIRTUAL_PARTICLE) {
	for (uint i = 0; i < cell_size; ++i) {
	    // particle number of cell placeholder
	    const int m = s_n[i];
	    // skip placeholder particles
	    if (m == VIRTUAL_PARTICLE)
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
template <uint cell_size>
__global__ void update_neighbours(int const* g_cell, int* g_nbl, float4* __empty__)
{
    // load particle from cell placeholder
    const int n = g_cell[GTID];
    const float3 r = unpack(tex1Dfetch(tex<float4>::r, n));
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

    // visit this cell
    update_cell_neighbours<cell_size, float4, true>(make_int3( 0,  0,  0), g_cell, g_nbl, r, n, count);
    // visit 26 neighbour cells, grouped into 13 pairs of mutually opposite cells
    update_cell_neighbours<cell_size, float4, false>(make_int3(-1, -1, -1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, float4, false>(make_int3(+1, +1, +1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, float4, false>(make_int3(-1, -1, +1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, float4, false>(make_int3(+1, +1, -1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, float4, false>(make_int3(-1, +1, +1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, float4, false>(make_int3(+1, -1, -1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, float4, false>(make_int3(+1, -1, +1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, float4, false>(make_int3(-1, +1, -1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, float4, false>(make_int3(-1, -1,  0), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, float4, false>(make_int3(+1, +1,  0), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, float4, false>(make_int3(-1, +1,  0), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, float4, false>(make_int3(+1, -1,  0), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, float4, false>(make_int3(-1,  0, -1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, float4, false>(make_int3(+1,  0, +1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, float4, false>(make_int3(-1,  0, +1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, float4, false>(make_int3(+1,  0, -1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, float4, false>(make_int3( 0, -1, -1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, float4, false>(make_int3( 0, +1, +1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, float4, false>(make_int3( 0, -1, +1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, float4, false>(make_int3( 0, +1, -1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, float4, false>(make_int3(-1,  0,  0), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, float4, false>(make_int3(+1,  0,  0), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, float4, false>(make_int3( 0, -1,  0), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, float4, false>(make_int3( 0, +1,  0), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, float4, false>(make_int3( 0,  0, -1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, float4, false>(make_int3( 0,  0, +1), g_cell, g_nbl, r, n, count);
}

template <uint cell_size>
__global__ void update_neighbours(int const* g_cell, int* g_nbl, float2* __empty__)
{
    // load particle from cell placeholder
    const int n = g_cell[GTID];
    const float2 r = unpack(tex1Dfetch(tex<float2>::r, n));
    // number of particles in neighbour list
    uint count = 0;

    // visit this cell
    update_cell_neighbours<cell_size, float2, true>(make_int2( 0,  0), g_cell, g_nbl, r, n, count);
    // visit 8 neighbour cells, grouped into 4 pairs of mutually opposite cells
    update_cell_neighbours<cell_size, float2, false>(make_int2(-1, -1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, float2, false>(make_int2(+1, +1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, float2, false>(make_int2(-1, +1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, float2, false>(make_int2(+1, -1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, float2, false>(make_int2(-1,  0), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, float2, false>(make_int2(+1,  0), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, float2, false>(make_int2( 0, -1), g_cell, g_nbl, r, n, count);
    update_cell_neighbours<cell_size, float2, false>(make_int2( 0, +1), g_cell, g_nbl, r, n, count);
}

/**
 * compute cell indices for given particle positions
 */
static __global__ void compute_cell(float4 const* g_part, uint* g_cell)
{
    float3 r = unpack(g_part[GTID]);
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
    g_cell[GTID] = uint(r.x) + ncell * (uint(r.y) + ncell * uint(r.z));
}

static __global__ void compute_cell(float2 const* g_part, uint* g_cell)
{
    float2 r = unpack(g_part[GTID]);
    r = (__saturatef(r / box) * (1.f - FLT_EPSILON)) * ncell;
    g_cell[GTID] = uint(r.x) + ncell * uint(r.y);
}

/**
 * compute global cell offsets in particle list
 */
static __global__ void find_cell_offset(uint* g_cell, int* g_cell_offset)
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
    g_otag[blockIdx.x * cell_size + threadIdx.x] = tag;
}

/**
 * generate ascending index sequence
 */
static __global__ void gen_index(int* g_idx)
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
    g_or[GTID] = tex1Dfetch(tex<U>::r, j);
    g_oR[GTID] = tex1Dfetch(tex<U>::R, j);
    g_ov[GTID] = tex1Dfetch(tex<U>::v, j);
    // permute particle tracking number
    g_otag[GTID] = tex1Dfetch(tex<>::tag, j);
}

/**
 * device constant wrappers
 */
cuda::symbol<uint> ljfluid_base<ljfluid_impl_gpu_neighbour>::npart(gpu::npart);
cuda::symbol<float> ljfluid_base<ljfluid_impl_gpu_neighbour>::box(gpu::box);
cuda::symbol<float> ljfluid_base<ljfluid_impl_gpu_neighbour>::timestep(gpu::timestep);
cuda::symbol<float> ljfluid_base<ljfluid_impl_gpu_neighbour>::r_cut(gpu::r_cut);
cuda::symbol<float> ljfluid_base<ljfluid_impl_gpu_neighbour>::rr_cut(gpu::rr_cut);
cuda::symbol<float> ljfluid_base<ljfluid_impl_gpu_neighbour>::en_cut(gpu::en_cut);
cuda::symbol<float> ljfluid_base<ljfluid_impl_gpu_neighbour>::rri_smooth(gpu::rri_smooth);

cuda::symbol<uint> ljfluid_base<ljfluid_impl_gpu_neighbour>::ncell(gpu::ncell);
cuda::symbol<uint> ljfluid_base<ljfluid_impl_gpu_neighbour>::nbl_size(gpu::nbl_size);
cuda::symbol<uint> ljfluid_base<ljfluid_impl_gpu_neighbour>::nbl_stride(gpu::nbl_stride);
cuda::symbol<float> ljfluid_base<ljfluid_impl_gpu_neighbour>::r_cell(gpu::r_cell);
cuda::symbol<float> ljfluid_base<ljfluid_impl_gpu_neighbour>::rr_cell(gpu::rr_cell);

/**
 * device texture wrappers
 */
cuda::texture<int> ljfluid_base<ljfluid_impl_gpu_neighbour>::tag(gpu::tex<>::tag);
cuda::texture<float4> ljfluid<ljfluid_impl_gpu_neighbour<3> >::r(gpu::tex<float4>::r);
cuda::texture<float4> ljfluid<ljfluid_impl_gpu_neighbour<3> >::R(gpu::tex<float4>::R);
cuda::texture<float4> ljfluid<ljfluid_impl_gpu_neighbour<3> >::v(gpu::tex<float4>::v);
cuda::texture<float2> ljfluid<ljfluid_impl_gpu_neighbour<2> >::r(gpu::tex<float2>::r);
cuda::texture<float2> ljfluid<ljfluid_impl_gpu_neighbour<2> >::R(gpu::tex<float2>::R);
cuda::texture<float2> ljfluid<ljfluid_impl_gpu_neighbour<2> >::v(gpu::tex<float2>::v);

/**
 * device function wrappers
 */
cuda::function<void (float3*, const float2)>
    ljfluid_base<ljfluid_impl_gpu_neighbour>::sample_smooth_function(gpu::sample_smooth_function);
cuda::function<void (int*)>
    ljfluid_base<ljfluid_impl_gpu_neighbour>::init_tags(gpu::init_tags);
cuda::function<void (uint const*, int const*, int const*, int*)>
    ljfluid_base<ljfluid_impl_gpu_neighbour>::assign_cells(gpu::assign_cells<CELL_SIZE>);
cuda::function<void (uint*, int*)>
    ljfluid_base<ljfluid_impl_gpu_neighbour>::find_cell_offset(gpu::find_cell_offset);
cuda::function<void (int*)>
    ljfluid_base<ljfluid_impl_gpu_neighbour>::gen_index(gpu::gen_index);

cuda::function<void (float4*, float4*, float4*, float4 const*)>
    ljfluid<ljfluid_impl_gpu_neighbour<3> >::inteq(gpu::inteq<float3>);
cuda::function<void (float4 const*, float4*, float4*, int const*, float*, float*)>
    ljfluid<ljfluid_impl_gpu_neighbour<3> >::mdstep(gpu::mdstep<float3, dfloat3>);
cuda::function<void (int const*, int*, float4*)>
    ljfluid<ljfluid_impl_gpu_neighbour<3> >::update_neighbours(gpu::update_neighbours<CELL_SIZE>);
cuda::function<void (float4 const*, uint*)>
    ljfluid<ljfluid_impl_gpu_neighbour<3> >::compute_cell(gpu::compute_cell);
cuda::function<void (const int*, float4*, float4*, float4*, int*)>
    ljfluid<ljfluid_impl_gpu_neighbour<3> >::order_particles(gpu::order_particles);

cuda::function<void (float2*, float2*, float2*, float2 const*)>
    ljfluid<ljfluid_impl_gpu_neighbour<2> >::inteq(gpu::inteq<float2>);
cuda::function<void (float2 const*, float2*, float2*, int const*, float*, float*)>
    ljfluid<ljfluid_impl_gpu_neighbour<2> >::mdstep(gpu::mdstep<float2, dfloat2>);
cuda::function<void (int const*, int*, float2*)>
    ljfluid<ljfluid_impl_gpu_neighbour<2> >::update_neighbours(gpu::update_neighbours<CELL_SIZE>);
cuda::function<void (float2 const*, uint*)>
    ljfluid<ljfluid_impl_gpu_neighbour<2> >::compute_cell(gpu::compute_cell);
cuda::function<void (const int*, float2*, float2*, float2*, int*)>
    ljfluid<ljfluid_impl_gpu_neighbour<2> >::order_particles(gpu::order_particles);

}} // namespace ljgpu::gpu
