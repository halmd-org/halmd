/*
 * Copyright © 2008-2011  Peter Colberg
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

#include <halmd/mdsim/gpu/box_kernel.cuh>
#include <halmd/mdsim/gpu/neighbours/from_binning_kernel.hpp>
#include <halmd/mdsim/gpu/particle_kernel.cuh>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/utility/gpu/thread.cuh>
#include <halmd/utility/gpu/variant.cuh>

using namespace halmd::mdsim::gpu::particle_kernel;
using namespace halmd::utility::gpu;

namespace halmd {
namespace mdsim {
namespace gpu {
namespace neighbours {
namespace from_binning_kernel {

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
/** positions, tags */
texture<float4> r_;

/**
 * compute neighbour cell
 */
__device__ unsigned int compute_neighbour_cell(fixed_vector<int, 3> const &offset)
{
    fixed_vector<int, 3> ncell(static_cast<fixed_vector<unsigned int, 3> >(get<3>(ncell_)));
    fixed_vector<int, 3> cell;

    // cell belonging to this execution block
    cell[0] = BID % ncell[0];
    cell[1] = (BID / ncell[0]) % ncell[1];
    cell[2] = BID / ncell[0] / ncell[1];
    // neighbour cell of this cell
    cell = element_mod(cell + ncell + offset, ncell);

    return (cell[2] * ncell[1] + cell[1]) * ncell[0] + cell[0];
}

__device__ unsigned int compute_neighbour_cell(fixed_vector<int, 2> const& offset)
{
    fixed_vector<int, 2> ncell(static_cast<fixed_vector<unsigned int, 2> >(get<2>(ncell_)));
    fixed_vector<int, 2> cell;

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
    tie(s_r[threadIdx.x], s_type[threadIdx.x]) = untagged<fixed_vector<float, dimension> >(tex1Dfetch(r_, n_));
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
    fixed_vector<float, dimension> r;
    tie(r, type) = untagged<fixed_vector<float, dimension> >(tex1Dfetch(r_, n));
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

    fixed_vector<int, dimension> j;
    for (j[0] = -1; j[0] <= 1; ++j[0]) {
        for (j[1] = -1; j[1] <= 1; ++j[1]) {
            if (dimension == 3) {
                for (j[2] = -1; j[2] <= 1; ++j[2]) {
                    if (j[0] == 0 && j[1] == 0 && j[2] == 0) {
                        goto self;
                    }
                    // visit 26 neighbour cells, grouped into 13 pairs of mutually opposite cells
                    update_cell_neighbours<false>(j, g_cell, r, type, n, count, g_neighbour);
                    update_cell_neighbours<false>(-j, g_cell, r, type, n, count, g_neighbour);
                }
            }
            else {
                if (j[0] == 0 && j[1] == 0) {
                    goto self;
                }
                // visit 8 neighbour cells, grouped into 4 pairs of mutually opposite cells
                update_cell_neighbours<false>(j, g_cell, r, type, n, count, g_neighbour);
                update_cell_neighbours<false>(-j, g_cell, r, type, n, count, g_neighbour);
            }
        }
    }

self:
    update_cell_neighbours<true>(j, g_cell, r, type, n, count, g_neighbour);

    // return failure if any neighbour list is fully occupied
    if (count == neighbour_size_) {
        *g_ret = EXIT_FAILURE;
    }
}

} // namespace from_binning_kernel

template <int dimension>
from_binning_wrapper<dimension> from_binning_wrapper<dimension>::kernel = {
    from_binning_kernel::rr_cut_skin_
  , get<dimension>(from_binning_kernel::ncell_)
  , from_binning_kernel::neighbour_size_
  , from_binning_kernel::neighbour_stride_
  , from_binning_kernel::nbox_
  , from_binning_kernel::r_
  , get<dimension>(from_binning_kernel::box_length_)
  , from_binning_kernel::update_neighbours<dimension>
};

template class from_binning_wrapper<3>;
template class from_binning_wrapper<2>;

} // namespace neighbours
} // namespace gpu
} // namespace mdsim
} // namespace halmd
