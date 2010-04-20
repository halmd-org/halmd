/*
 * Copyright © 2008-2010  Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
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

#include <boost/mpl/if.hpp>
#include <float.h>

#include <halmd/algorithm/gpu/base.cuh>
#include <halmd/mdsim/gpu/box_kernel.cuh>
#include <halmd/mdsim/gpu/particle_kernel.cuh>
#include <halmd/numeric/gpu/blas/symmetric.cuh>
#include <halmd/numeric/gpu/blas/vector.cuh>

using namespace boost::mpl;
using namespace halmd::mdsim::gpu::particle_kernel;
using namespace halmd::numeric::gpu::blas;

namespace halmd { namespace mdsim { namespace gpu
{

namespace neighbor_kernel
{

/** (cutoff lengths + neighbor list skin)² */
texture<float, 1, cudaReadModeElementType> rr_cut_skin_;
/** neighbor list length */
__constant__ unsigned int neighbor_size_;
/** neighbor list stride */
__constant__ unsigned int neighbor_stride_;
/** number of particles in simulation box */
__constant__ unsigned int nbox_;

template <size_t N>
struct dim_
{
    /** positions, tags */
    static texture<float4, 1, cudaReadModeElementType> r;
    /** number of cells per dimension */
    static __constant__ typename if_c<N == 3, uint3, uint2>::type ncell;
    /** cubic box edgle length */
    static __constant__ typename if_c<N == 3, float3, float2>::type box_length;
    /** cell edge lengths */
    static __constant__ typename if_c<N == 3, float3, float2>::type cell_length;
};

// explicit instantiation
template class dim_<3>;
template class dim_<2>;

/**
 * compute neighbor cell
 */
__device__ unsigned int compute_neighbor_cell(int3 const &offset)
{
    vector<unsigned int, 3> ncell = dim_<3>::ncell;
    // cell belonging to this execution block
    int x = BID % ncell[0];
    int y = (BID / ncell[0]) % ncell[1];
    int z = BID / ncell[0] / ncell[1];
    // neighbor cell of this cell
    x = (x + ncell[0] + offset.x) % ncell[0];
    y = (y + ncell[1] + offset.y) % ncell[1];
    z = (z + ncell[2] + offset.z) % ncell[2];

    return (z * ncell[1] + y) * ncell[0] + x;
}

__device__ unsigned int compute_neighbor_cell(int2 const& offset)
{
    vector<unsigned int, 2> ncell = dim_<2>::ncell;
    // cell belonging to this execution block
    int x = BID % ncell[0];
    int y = BID / ncell[0];
    // neighbor cell of this cell
    x = (x + ncell[0] + offset.x) % ncell[0];
    y = (y + ncell[1] + offset.y) % ncell[1];

    return y * ncell[0] + x;
}

/**
 * update neighbor list with particles of given cell
 */
template <bool same_cell, typename T, typename I>
__device__ void update_cell_neighbors(
  I const& offset,
  unsigned int const* g_cell,
  T const& r,
  unsigned int type,
  unsigned int const& n,
  unsigned int& count,
  unsigned int* g_neighbor)
{
    extern __shared__ unsigned int s_n[];
    unsigned int* const s_type = &s_n[blockDim.x];
    T* const s_r = reinterpret_cast<T*>(&s_n[2 * blockDim.x]);
    enum { dimension = T::static_size };

    // shared memory barrier
    __syncthreads();

    // compute cell index
    unsigned int const cell = compute_neighbor_cell(offset);
    // load particles in cell
    unsigned int const n_ = g_cell[cell * blockDim.x + threadIdx.x];
    s_n[threadIdx.x] = n_;
    s_r[threadIdx.x] = untagged<vector<float, dimension> >(tex1Dfetch(dim_<dimension>::r, n_), s_type[threadIdx.x]);
    __syncthreads();

    if (n == PLACEHOLDER) return;

    for (unsigned int i = 0; i < blockDim.x; ++i) {
        // particle number of cell placeholder
        unsigned int const m = s_n[i];
        // skip placeholder particles
        if (m == PLACEHOLDER) break;
        // skip same particle
        if (same_cell && i == threadIdx.x) continue;

        // particle distance vector
        T dr = r - s_r[i];
        // enforce periodic boundary conditions
        T L = dim_<dimension>::box_length;
        box_kernel::reduce_periodic(dr, L);
        // squared particle distance
        float rr = inner_prod(dr, dr);

        // enforce cutoff length with neighbor list skin
        float rr_cut_skin = tex1Dfetch(rr_cut_skin_, symmetric_matrix::lower_index(type, s_type[i]));
        if (rr <= rr_cut_skin && count < neighbor_size_) {
            // scattered write to neighbor list
            g_neighbor[count * neighbor_stride_ + n] = m;
            // increment neighbor list particle count
            count++;
        }
    }
}

/**
 * update neighbor lists
 */
template <unsigned int dimension>
__global__ void update_neighbors(
  unsigned int* g_neighbor,
  unsigned int* g_ret,
  unsigned int const* g_cell)
{
    // load particle from cell placeholder
    unsigned int const n = g_cell[GTID];
    unsigned int type;
    vector<float, dimension> const r = untagged<vector<float, dimension> >(tex1Dfetch(dim_<dimension>::r, n), type);
    // number of particles in neighbor list
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
    // single neighbor cell will more or less have the same direction.
    // Thus, when adding the force sums of all neighbor cells, we add
    // huge force sums which will mostly cancel each other out in an
    // equilibrated system, giving a small and very inaccurate total
    // force due to being limited to single-precision floating-point
    // arithmetic.
    //
    // Besides implementing the summation in double precision arithmetic,
    // choosing the order of summation over cells such that one partial
    // neighbor cell force sum is always followed by the sum of the
    // opposite neighbor cell softens the velocity drift.
    //

    if (dimension == 3) {
        // visit this cell
        update_cell_neighbors<true>(make_int3( 0,  0,  0), g_cell, r, type, n, count, g_neighbor);
        // visit 26 neighbor cells, grouped into 13 pairs of mutually opposite cells
        update_cell_neighbors<false>(make_int3(-1, -1, -1), g_cell, r, type, n, count, g_neighbor);
        update_cell_neighbors<false>(make_int3(+1, +1, +1), g_cell, r, type, n, count, g_neighbor);
        update_cell_neighbors<false>(make_int3(-1, -1, +1), g_cell, r, type, n, count, g_neighbor);
        update_cell_neighbors<false>(make_int3(+1, +1, -1), g_cell, r, type, n, count, g_neighbor);
        update_cell_neighbors<false>(make_int3(-1, +1, +1), g_cell, r, type, n, count, g_neighbor);
        update_cell_neighbors<false>(make_int3(+1, -1, -1), g_cell, r, type, n, count, g_neighbor);
        update_cell_neighbors<false>(make_int3(+1, -1, +1), g_cell, r, type, n, count, g_neighbor);
        update_cell_neighbors<false>(make_int3(-1, +1, -1), g_cell, r, type, n, count, g_neighbor);
        update_cell_neighbors<false>(make_int3(-1, -1,  0), g_cell, r, type, n, count, g_neighbor);
        update_cell_neighbors<false>(make_int3(+1, +1,  0), g_cell, r, type, n, count, g_neighbor);
        update_cell_neighbors<false>(make_int3(-1, +1,  0), g_cell, r, type, n, count, g_neighbor);
        update_cell_neighbors<false>(make_int3(+1, -1,  0), g_cell, r, type, n, count, g_neighbor);
        update_cell_neighbors<false>(make_int3(-1,  0, -1), g_cell, r, type, n, count, g_neighbor);
        update_cell_neighbors<false>(make_int3(+1,  0, +1), g_cell, r, type, n, count, g_neighbor);
        update_cell_neighbors<false>(make_int3(-1,  0, +1), g_cell, r, type, n, count, g_neighbor);
        update_cell_neighbors<false>(make_int3(+1,  0, -1), g_cell, r, type, n, count, g_neighbor);
        update_cell_neighbors<false>(make_int3( 0, -1, -1), g_cell, r, type, n, count, g_neighbor);
        update_cell_neighbors<false>(make_int3( 0, +1, +1), g_cell, r, type, n, count, g_neighbor);
        update_cell_neighbors<false>(make_int3( 0, -1, +1), g_cell, r, type, n, count, g_neighbor);
        update_cell_neighbors<false>(make_int3( 0, +1, -1), g_cell, r, type, n, count, g_neighbor);
        update_cell_neighbors<false>(make_int3(-1,  0,  0), g_cell, r, type, n, count, g_neighbor);
        update_cell_neighbors<false>(make_int3(+1,  0,  0), g_cell, r, type, n, count, g_neighbor);
        update_cell_neighbors<false>(make_int3( 0, -1,  0), g_cell, r, type, n, count, g_neighbor);
        update_cell_neighbors<false>(make_int3( 0, +1,  0), g_cell, r, type, n, count, g_neighbor);
        update_cell_neighbors<false>(make_int3( 0,  0, -1), g_cell, r, type, n, count, g_neighbor);
        update_cell_neighbors<false>(make_int3( 0,  0, +1), g_cell, r, type, n, count, g_neighbor);
    }
    else {
        // visit this cell
        update_cell_neighbors<true>(make_int2( 0,  0), g_cell, r, type, n, count, g_neighbor);
        // visit 8 neighbor cells, grouped into 4 pairs of mutually opposite cells
        update_cell_neighbors<false>(make_int2(-1, -1), g_cell, r, type, n, count, g_neighbor);
        update_cell_neighbors<false>(make_int2(+1, +1), g_cell, r, type, n, count, g_neighbor);
        update_cell_neighbors<false>(make_int2(-1, +1), g_cell, r, type, n, count, g_neighbor);
        update_cell_neighbors<false>(make_int2(+1, -1), g_cell, r, type, n, count, g_neighbor);
        update_cell_neighbors<false>(make_int2(-1,  0), g_cell, r, type, n, count, g_neighbor);
        update_cell_neighbors<false>(make_int2(+1,  0), g_cell, r, type, n, count, g_neighbor);
        update_cell_neighbors<false>(make_int2( 0, -1), g_cell, r, type, n, count, g_neighbor);
        update_cell_neighbors<false>(make_int2( 0, +1), g_cell, r, type, n, count, g_neighbor);
    }

    // return failure if any neighbor list is fully occupied
    if (count == neighbor_size_) {
        *g_ret = EXIT_FAILURE;
    }
}

/**
 * compute cell indices for given particle positions
 */
template <typename vector_type>
__device__ inline unsigned int compute_cell_index(vector_type r)
{
    enum { dimension = vector_type::static_size };
    vector_type L = dim_<dimension>::box_length;
    vector<unsigned int, dimension> ncell = dim_<dimension>::ncell;

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
    vector_type frac = __saturate(element_div(r, L)) * (1.f - FLT_EPSILON);
    vector<unsigned int, dimension> index(element_prod(static_cast<vector_type>(ncell), frac));
    if (dimension == 3) {
        return index[0] + ncell[0] * (index[1] + ncell[1] * index[2]);
    }
    else {
        return index[0] + ncell[0] * index[1];
    }
}

/**
 * compute cell indices for particle positions
 */
template <unsigned int dimension>
__global__ void compute_cell(float4 const* g_r, unsigned int* g_cell)
{
    vector<float, dimension> const r = untagged<vector<float, dimension> >(g_r[GTID]);
    g_cell[GTID] = compute_cell_index(r);
}

/**
 * compute global cell offsets in particle list
 */
__global__ void find_cell_offset(unsigned int* g_cell, unsigned int* g_cell_offset)
{
    const unsigned int j = g_cell[GTID];
    const unsigned int k = (GTID > 0 && GTID < nbox_) ? g_cell[GTID - 1] : j;

    if (GTID == 0 || k < j) {
        // particle marks the start of a cell
        g_cell_offset[j] = GTID;
    }
}

/**
 * assign particles to cells
 */
__global__ void assign_cells(
  unsigned int* g_ret,
  unsigned int const* g_cell,
  unsigned int const* g_cell_offset,
  unsigned int const* g_itag,
  unsigned int* g_otag)
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
    unsigned int tag = PLACEHOLDER;
    // mark as real particle if appropriate
    if (offset != PLACEHOLDER && n < nbox_ && g_cell[n] == BID) {
        tag = g_itag[n];
    }
    // return failure if any cell list is fully occupied
    if (tag != PLACEHOLDER && (threadIdx.x + 1) == blockDim.x) {
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
    g_index[GTID] = (GTID < nbox_) ? GTID : 0;
}

} // namespace neighbor_kernel

}}} //namespace halmd::mdsim::gpu
