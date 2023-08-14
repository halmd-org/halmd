/*
 * Copyright © 2021      Jaslo Ziska
 * Copyright © 2008-2011 Peter Colberg
 * Copyright © 2014      Nicolas Höft
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <halmd/mdsim/gpu/box_kernel.cuh>
#include <halmd/mdsim/gpu/neighbours/from_binning_kernel.hpp>
#include <halmd/mdsim/gpu/particle_kernel.cuh>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/utility/gpu/thread.cuh>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace neighbours {
namespace from_binning_kernel {

/**
 * compute neighbour cell
 */
inline __device__ unsigned int compute_neighbour_cell(
    unsigned int const cell_index
  , fixed_vector<int, 3> const& offset
  , fixed_vector<int, 3> const& ncell
)
{
    // cell index of specified cell offset
    fixed_vector<int, 3> cell;
    cell[0] = cell_index % ncell[0];
    cell[1] = (cell_index / ncell[0]) % ncell[1];
    cell[2] = cell_index / ncell[0] / ncell[1];
    // neighbour cell of this cell
    cell = element_mod(cell + ncell + offset, ncell);

    return (cell[2] * ncell[1] + cell[1]) * ncell[0] + cell[0];
}

inline __device__ unsigned int compute_neighbour_cell(
    unsigned int const cell_index
  , fixed_vector<int, 2> const& offset
  , fixed_vector<int, 2> const& ncell
)
{
    // cell index of specified cell offset
    fixed_vector<int, 2> cell;
    cell[0] = cell_index % ncell[0];
    cell[1] = cell_index / ncell[0];
    // neighbour cell of this cell
    cell = element_mod(cell + ncell + offset, ncell);

    return cell[1] * ncell[0] + cell[0];
}

/**
 * compute cell indices for given particle positions
 */
template <typename vector_type, typename cell_size_type>
inline __device__ unsigned int compute_cell_index(
    vector_type r
  , vector_type cell_length
  , cell_size_type ncell
)
{
    enum { dimension = vector_type::static_size };

    cell_size_type index = element_mod(
        static_cast<cell_size_type>(element_div(r, cell_length) + static_cast<vector_type>(ncell))
      , ncell
    );
    // FIXME check PTX to ensure CUDA unrolls this loop
    unsigned int offset = index[dimension - 1];
    for (int i = dimension - 2; i >= 0; i--) {
        offset *= ncell[i];
        offset += index[i];
    }
    return offset;
}

/**
 * update neighbour list with particles of given cell
 */
template <bool same_cell, bool unroll_force_loop, typename vector_type, typename cell_size_type, typename cell_difference_type>
__device__ void update_cell_neighbours(
    cudaTextureObject_t t_rr_cut_skin
  , cudaTextureObject_t t_r2
  , cell_difference_type const& offset
  , cell_size_type const& ncell
  , unsigned int const* g_cell1
  , unsigned int const* g_cell2
  , vector_type const& r
  , unsigned int type
  , unsigned int ntype1
  , unsigned int ntype2
  , unsigned int const& n
  , unsigned int& count
  , unsigned int* g_neighbour
  , unsigned int neighbour_size
  , unsigned int neighbour_stride
  , vector_type const& box_length
)
{
    extern __shared__ unsigned int s_n[];
    unsigned int* const s_type = &s_n[blockDim.x];
    vector_type* const s_r = reinterpret_cast<vector_type*>(&s_n[2 * blockDim.x]);

    // shared memory barrier
    __syncthreads();

    // compute cell index
    unsigned int const cell = compute_neighbour_cell(BID, offset, static_cast<cell_difference_type>(ncell));
    // load particles in cell
    unsigned int const n_ = g_cell2[cell * blockDim.x + threadIdx.x];
    s_n[threadIdx.x] = n_;
    tie(s_r[threadIdx.x], s_type[threadIdx.x]) <<= tex1Dfetch<float4>(t_r2, n_);
    __syncthreads();

    if (n == particle_kernel::placeholder) return;

    for (unsigned int i = 0; i < blockDim.x; ++i) {
        // particle number of cell placeholder
        unsigned int const m = s_n[i];
        // skip placeholder particles
        if (m == particle_kernel::placeholder) break;
        // skip same particle
        if (same_cell && n == m && g_cell1 == g_cell2) continue;

        // particle distance vector
        vector_type dr = r - s_r[i];
        // enforce periodic boundary conditions
        box_kernel::reduce_periodic(dr, box_length);
        // squared particle distance
        float rr = inner_prod(dr, dr);

        // enforce cutoff distance with neighbour list skin
        float rr_cut_skin = tex1Dfetch<float>(t_rr_cut_skin, type * ntype2 + s_type[i]);

        if (rr <= rr_cut_skin && count < neighbour_size) {
            // scattered write to neighbour list
            if (unroll_force_loop) {
                g_neighbour[n * neighbour_size + count] = m;
            } else {
                g_neighbour[count * neighbour_stride + n] = m;
            }
            // increment neighbour list particle count
            count++;
        }
    }
}

/**
 * update neighbour lists
 */
template <bool unroll_force_loop, unsigned int dimension>
__global__ void update_neighbours(
    cudaTextureObject_t t_rr_cut_skin
  , cudaTextureObject_t t_r1
  , cudaTextureObject_t t_r2
  , int* g_ret
  , unsigned int* g_neighbour
  , unsigned int neighbour_size
  , unsigned int neighbour_stride
  , unsigned int const* g_cell1
  , unsigned int const* g_cell2
  , unsigned int ntype1
  , unsigned int ntype2
  , unsigned int total_cell_size1
  , fixed_vector<unsigned int, dimension> ncell
  , fixed_vector<float, dimension> box_length
)
{
    // load particle from cell placeholder
    unsigned int const n = g_cell1[GTID];
    unsigned int type;
    fixed_vector<float, dimension> r;
    tie(r, type) <<= tex1Dfetch<float4>(t_r1, n);
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
    // algorithm, the centre of mass velocity effectively stays zero if
    // the initial list of particles arranged on a lattice is randomly
    // permuted before simulation.
    // Using the cell algorithm with the WCA potential however results
    // in a continuously drifting centre of mass velocity, independent
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

    fixed_vector<int, dimension> j;
    for (j[0] = -1; j[0] <= 1; ++j[0]) {
        for (j[1] = -1; j[1] <= 1; ++j[1]) {
            if (dimension == 3) {
                for (j[2] = -1; j[2] <= 1; ++j[2]) {
                    if (j[0] == 0 && j[1] == 0 && j[2] == 0) {
                        goto self;
                    }
                    // visit 26 neighbour cells, grouped into 13 pairs of mutually opposite cells
                    update_cell_neighbours<false, unroll_force_loop>(t_rr_cut_skin, t_r2, j, ncell, g_cell1, g_cell2, r, type, ntype1,
                        ntype2, n, count, g_neighbour, neighbour_size, neighbour_stride, box_length);
                    update_cell_neighbours<false, unroll_force_loop>(t_rr_cut_skin, t_r2, -j, ncell, g_cell1, g_cell2, r, type, ntype1,
                        ntype2, n, count, g_neighbour, neighbour_size, neighbour_stride, box_length);
                }
            }
            else {
                if (j[0] == 0 && j[1] == 0) {
                    goto self;
                }
                // visit 8 neighbour cells, grouped into 4 pairs of mutually opposite cells
                update_cell_neighbours<false, unroll_force_loop>(t_rr_cut_skin, t_r2, j, ncell, g_cell1, g_cell2, r, type, ntype1,
                    ntype2, n, count, g_neighbour, neighbour_size, neighbour_stride, box_length);

                update_cell_neighbours<false, unroll_force_loop>(t_rr_cut_skin, t_r2, -j, ncell, g_cell1, g_cell2, r, type, ntype1,
                    ntype2, n, count, g_neighbour, neighbour_size, neighbour_stride, box_length);
            }
        }
    }

self:
    update_cell_neighbours<true, unroll_force_loop>(t_rr_cut_skin, t_r2, j, ncell, g_cell1, g_cell2, r, type, ntype1,
        ntype2, n, count, g_neighbour, neighbour_size, neighbour_stride, box_length);

    // return failure if any neighbour list is fully occupied
    if (count == neighbour_size) {
        *g_ret = EXIT_FAILURE;
    }
}


/**
 * update neighbour list with particles of given cell using a naive implementation
 * that iterates through the given neighbour cell
 * As the threads do not necessarily have particles that are in the same cell, we hope
 * that the cache hides this bottleneck.
 */
template <bool same_cell, bool unroll_force_loop, typename vector_type, typename cell_size_type, typename cell_difference_type>
__device__ void update_cell_neighbours_naive(
    cudaTextureObject_t t_rr_cut_skin
  , cudaTextureObject_t t_r2
  ,  cell_difference_type const& offset
  , unsigned int const cell_index
  , cell_size_type const& ncell
  , unsigned int const* g_cell
  , bool same_particle
  , vector_type const& r
  , unsigned int type
  , unsigned int ntype1
  , unsigned int ntype2
  , unsigned int const& n
  , unsigned int& count
  , unsigned int* g_neighbour
  , unsigned int neighbour_size
  , unsigned int neighbour_stride
  , unsigned int cell_size
  , vector_type const& box_length
)
{
    // compute cell index
    unsigned int const cell = compute_neighbour_cell(cell_index, offset, static_cast<cell_difference_type>(ncell));

    for (unsigned int i = 0; i < cell_size; ++i) {
        // load index of particle in neighbour cell
        unsigned int const m = g_cell[cell * cell_size + i];
        // skip placeholder particles
        if (m == particle_kernel::placeholder) break;

        // skip same particle
        if (same_cell && m == n && same_particle) continue;

        vector_type r2;
        unsigned int type2;
        tie(r2, type2) <<= tex1Dfetch<float4>(t_r2, m);

        // particle distance vector
        vector_type dr = r - r2;
        // enforce periodic boundary conditions
        box_kernel::reduce_periodic(dr, box_length);
        // squared particle distance
        float rr = inner_prod(dr, dr);

        // enforce cutoff distance with neighbour list skin
        float rr_cut_skin = tex1Dfetch<float>(t_rr_cut_skin, type * ntype2 + type2);

        if (rr <= rr_cut_skin && count < neighbour_size) {
            // scattered write to neighbour list
            if (unroll_force_loop) {
                g_neighbour[n * neighbour_size + count] = m;
            } else {
                g_neighbour[count * neighbour_stride + n] = m;
            }
            // increment neighbour list particle count
            count++;
        }
    }
}

/**
 * update neighbour lists
 *
 * This "naive" implementation calculates the cell index of particle "A" (whose
 * neighbour list is being constructed) on-the-fly and then iterates through the
 * neighbour cells of particle "B".
 * This can cause multiple independent reads of the same cell when A-particles in the
 * same block are not in the same cell.
 */
template <bool unroll_force_loop, unsigned int dimension>
__global__ void update_neighbours_naive(
    cudaTextureObject_t t_rr_cut_skin
  , cudaTextureObject_t t_r2
  , int* g_ret
  , float4 const* g_r1
  , unsigned int nparticle
  , bool same_particle
  , unsigned int* g_neighbour
  , unsigned int neighbour_size
  , unsigned int neighbour_stride
  , unsigned int const* g_cell2
  , unsigned int ntype1
  , unsigned int ntype2
  , fixed_vector<unsigned int, dimension> ncell
  , fixed_vector<float, dimension> cell_length
  , unsigned int cell_size
  , fixed_vector<float, dimension> box_length
)
{
    // make sure we do not read the position of particle placeholders
    if (GTID >= nparticle)
        return;
    unsigned int const n = GTID;
    // load particle from global memory associated with this thread
    unsigned int type;
    fixed_vector<float, dimension> r;
    tie(r, type) <<= g_r1[n];

    // number of particles in neighbour list
    unsigned int count = 0;
    // cell offset of particle
    unsigned int cell_index = compute_cell_index(r, cell_length, ncell);

    fixed_vector<int, dimension> j;
    for (j[0] = -1; j[0] <= 1; ++j[0]) {
        for (j[1] = -1; j[1] <= 1; ++j[1]) {
            if (dimension == 3) {
                for (j[2] = -1; j[2] <= 1; ++j[2]) {
                    if (j[0] == 0 && j[1] == 0 && j[2] == 0) {
                        goto self;
                    }
                    // visit 26 neighbour cells, grouped into 13 pairs of mutually opposite cells
                    update_cell_neighbours_naive<false, unroll_force_loop>(t_rr_cut_skin, t_r2, j, cell_index, ncell,
                        g_cell2, same_particle, r, type, ntype1, ntype2, n, count, g_neighbour, neighbour_size,
                        neighbour_stride, cell_size, box_length);
                    update_cell_neighbours_naive<false, unroll_force_loop>(t_rr_cut_skin, t_r2, -j, cell_index, ncell,
                        g_cell2, same_particle, r, type, ntype1, ntype2, n, count, g_neighbour, neighbour_size,
                        neighbour_stride, cell_size, box_length);
                }
            }
            else {
                if (j[0] == 0 && j[1] == 0) {
                    goto self;
                }
                // visit 8 neighbour cells, grouped into 4 pairs of mutually opposite cells
                update_cell_neighbours_naive<false, unroll_force_loop>(t_rr_cut_skin, t_r2, j, cell_index, ncell, g_cell2,
                    same_particle, r, type, ntype1, ntype2, n, count, g_neighbour, neighbour_size, neighbour_stride,
                    cell_size, box_length);
                update_cell_neighbours_naive<false, unroll_force_loop>(t_rr_cut_skin, t_r2, -j, cell_index, ncell,
                    g_cell2, same_particle, r, type, ntype1, ntype2, n, count, g_neighbour, neighbour_size,
                    neighbour_stride, cell_size, box_length);
            }
        }
    }

self:
    update_cell_neighbours_naive<true, unroll_force_loop>(t_rr_cut_skin, t_r2, j, cell_index, ncell, g_cell2,
        same_particle, r, type, ntype1, ntype2, n, count, g_neighbour, neighbour_size, neighbour_stride, cell_size,
        box_length);

    // return failure if any neighbour list is fully occupied
    if (count == neighbour_size) {
        *g_ret = EXIT_FAILURE;
    }
}

} // namespace from_binning_kernel

template <int dimension>
from_binning_wrapper<dimension> from_binning_wrapper<dimension>::kernel = {
    {
        from_binning_kernel::update_neighbours<true, dimension>
      , from_binning_kernel::update_neighbours_naive<true, dimension>
    }
  , {
        from_binning_kernel::update_neighbours<false, dimension>
      , from_binning_kernel::update_neighbours_naive<false, dimension>
    }
};

template class from_binning_wrapper<3>;
template class from_binning_wrapper<2>;

} // namespace neighbours
} // namespace gpu
} // namespace mdsim
} // namespace halmd
