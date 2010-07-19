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

#include <float.h>

#include <halmd/algorithm/gpu/reduction.cuh>
#include <halmd/mdsim/gpu/box_kernel.cuh>
#include <halmd/mdsim/gpu/neighbour_kernel.hpp>
#include <halmd/mdsim/gpu/particle_kernel.cuh>
#include <halmd/numeric/gpu/blas/symmetric.cuh>
#include <halmd/numeric/gpu/blas/vector.cuh>
#include <halmd/utility/gpu/thread.cuh>
#include <halmd/utility/gpu/variant.cuh>

using namespace halmd::algorithm::gpu;
using namespace halmd::mdsim::gpu::particle_kernel;
using namespace halmd::numeric::gpu::blas;
using namespace halmd::utility::gpu;

namespace halmd
{
namespace mdsim { namespace gpu
{
namespace neighbour_kernel
{

/** (cutoff lengths + neighbour list skin)² */
texture<float> rr_cut_skin_;
/** neighbour list length */
__constant__ unsigned int neighbour_size_;
/** neighbour list stride */
__constant__ unsigned int neighbour_stride_;
/** number of particles in simulation box */
__constant__ unsigned int nbox_;
/** cuboid box edgle length */
__constant__ variant<map<pair<int_<3>, float3>, pair<int_<2>, float2> > > box_length_;
/** number of cells per dimension */
__constant__ variant<map<pair<int_<3>, uint3>, pair<int_<2>, uint2> > > ncell_;
/** cell edge lengths */
__constant__ variant<map<pair<int_<3>, float3>, pair<int_<2>, float2> > > cell_length_;
/** positions, tags */
texture<float4> r_;

/**
 * compute neighbour cell
 */
__device__ unsigned int compute_neighbour_cell(vector<int, 3> const &offset)
{
    vector<int, 3> ncell(static_cast<vector<unsigned int, 3> >(get<3>(ncell_)));
    vector<int, 3> cell;

    // cell belonging to this execution block
    cell[0] = BID % ncell[0];
    cell[1] = (BID / ncell[0]) % ncell[1];
    cell[2] = BID / ncell[0] / ncell[1];
    // neighbour cell of this cell
    cell = element_mod(cell + ncell + offset, ncell);

    return (cell[2] * ncell[1] + cell[1]) * ncell[0] + cell[0];
}

__device__ unsigned int compute_neighbour_cell(vector<int, 2> const& offset)
{
    vector<int, 2> ncell(static_cast<vector<unsigned int, 2> >(get<2>(ncell_)));
    vector<int, 2> cell;

    // cell belonging to this execution block
    cell[0] = BID % ncell[0];
    cell[1] = BID / ncell[0];
    // neighbour cell of this cell
    cell = element_mod(cell + ncell + offset, ncell);

    return cell[1] * ncell[0] + cell[0];
}

/**
 * update neighbour list with particles of given cell
 */
template <bool same_cell, typename T, typename I>
__device__ void update_cell_neighbours(
  I const& offset,
  unsigned int const* g_cell,
  T const& r,
  unsigned int type,
  unsigned int const& n,
  unsigned int& count,
  unsigned int* g_neighbour)
{
    extern __shared__ unsigned int s_n[];
    unsigned int* const s_type = &s_n[blockDim.x];
    T* const s_r = reinterpret_cast<T*>(&s_n[2 * blockDim.x]);
    enum { dimension = T::static_size };

    // shared memory barrier
    __syncthreads();

    // compute cell index
    unsigned int const cell = compute_neighbour_cell(offset);
    // load particles in cell
    unsigned int const n_ = g_cell[cell * blockDim.x + threadIdx.x];
    s_n[threadIdx.x] = n_;
    tie(s_r[threadIdx.x], s_type[threadIdx.x]) = untagged<vector<float, dimension> >(tex1Dfetch(r_, n_));
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
        T L = get<dimension>(box_length_);
        box_kernel::reduce_periodic(dr, L);
        // squared particle distance
        float rr = inner_prod(dr, dr);

        // enforce cutoff length with neighbour list skin
        float rr_cut_skin = tex1Dfetch(rr_cut_skin_, symmetric_matrix::lower_index(type, s_type[i]));
        if (rr <= rr_cut_skin && count < neighbour_size_) {
            // scattered write to neighbour list
            g_neighbour[count * neighbour_stride_ + n] = m;
            // increment neighbour list particle count
            count++;
        }
    }
}

/**
 * update neighbour lists
 */
template <unsigned int dimension>
__global__ void update_neighbours(
  int* g_ret,
  unsigned int* g_neighbour,
  unsigned int const* g_cell)
{
    // load particle from cell placeholder
    unsigned int const n = g_cell[GTID];
    unsigned int type;
    vector<float, dimension> r;
    tie(r, type) = untagged<vector<float, dimension> >(tex1Dfetch(r_, n));
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

    for (int i = -1; i <= 1; ++i) {
        for (int j = -1; j <= 1; ++j) {
            if (dimension == 3) {
                for (int k = -1; k <= 1; ++k) {
                    if (i == 0 && j == 0 && k == 0) {
                        goto self;
                    }
                    // visit 26 neighbour cells, grouped into 13 pairs of mutually opposite cells
                    update_cell_neighbours<false>(make_int3(i, j, k), g_cell, r, type, n, count, g_neighbour);
                    update_cell_neighbours<false>(make_int3(-i, -j, -k), g_cell, r, type, n, count, g_neighbour);
                }
            }
            else {
                if (i == 0 && j == 0) {
                    goto self;
                }
                // visit 8 neighbour cells, grouped into 4 pairs of mutually opposite cells
                update_cell_neighbours<false>(make_int2(i, j), g_cell, r, type, n, count, g_neighbour);
                update_cell_neighbours<false>(make_int2(-i, -j), g_cell, r, type, n, count, g_neighbour);
            }
        }
    }

self:
    // visit this cell
    if (dimension == 3) {
        update_cell_neighbours<true>(make_int3( 0,  0,  0), g_cell, r, type, n, count, g_neighbour);
    }
    else {
        update_cell_neighbours<true>(make_int2( 0,  0), g_cell, r, type, n, count, g_neighbour);
    }

    // return failure if any neighbour list is fully occupied
    if (count == neighbour_size_) {
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
    vector_type L = get<dimension>(box_length_);
    vector<unsigned int, dimension> ncell = get<dimension>(ncell_);

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
    vector<float, dimension> r;
    unsigned int type;
    tie(r, type) = untagged<vector<float, dimension> >(g_r[GTID]);
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
  int* g_ret,
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

/**
 * maximum squared particle distance
 */
template <
    typename vector_type
  , int threads
>
__global__ void displacement(float4* g_r, float4* g_r0, typename vector_type::value_type* g_rr)
{
    typedef typename vector_type::value_type float_type;
    enum { dimension = vector_type::static_size };

    extern __shared__ char __s_array[]; // CUDA 3.0/3.1 breaks template __shared__ type
    float_type* const s_rr = reinterpret_cast<float_type*>(__s_array);
    float_type rr = 0;

    for (uint i = GTID; i < nbox_; i += GTDIM) {
        vector_type r;
        unsigned int type;
        tie(r, type) = untagged<vector_type>(g_r[i]);
        vector_type r0;
        tie(r0, type) = untagged<vector_type>(g_r0[i]);
        r -= r0;
        box_kernel::reduce_periodic(r, static_cast<vector_type>(get<dimension>(box_length_)));
        rr = max(rr, inner_prod(r, r));
    }

    // reduced values for this thread
    s_rr[TID] = rr;
    __syncthreads();

    // compute reduced value for all threads in block
    reduce<threads / 2, sum_>(rr, s_rr);

    if (TID < 1) {
        // store block reduced value in global memory
        g_rr[blockIdx.x] = rr;
    }
}

} // namespace neighbour_kernel

template <int dimension>
neighbour_wrapper<dimension> neighbour_wrapper<dimension>::kernel = {
    neighbour_kernel::rr_cut_skin_
  , get<dimension>(neighbour_kernel::ncell_)
  , neighbour_kernel::neighbour_size_
  , neighbour_kernel::neighbour_stride_
  , neighbour_kernel::nbox_
  , neighbour_kernel::r_
  , get<dimension>(neighbour_kernel::box_length_)
  , get<dimension>(neighbour_kernel::cell_length_)
  , neighbour_kernel::assign_cells
  , neighbour_kernel::find_cell_offset
  , neighbour_kernel::gen_index
  , neighbour_kernel::update_neighbours<dimension>
  , neighbour_kernel::compute_cell<dimension>
  , neighbour_kernel::displacement<vector<float, dimension>, 512>
  , neighbour_kernel::displacement<vector<float, dimension>, 256>
  , neighbour_kernel::displacement<vector<float, dimension>, 128>
  , neighbour_kernel::displacement<vector<float, dimension>, 64>
  , neighbour_kernel::displacement<vector<float, dimension>, 32>
};

template class neighbour_wrapper<3>;
template class neighbour_wrapper<2>;

}} //namespace mdsim::gpu

} //namespace halmd
