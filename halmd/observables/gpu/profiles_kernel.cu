/*
 * Copyright © 2010  Felix Höfling and Peter Colberg
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
#include <halmd/mdsim/force_kernel.hpp>
#include <halmd/mdsim/gpu/particle_kernel.cuh>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/observables/gpu/profiles_kernel.hpp>
#include <halmd/utility/gpu/thread.cuh>

using namespace halmd::algorithm::gpu;
using namespace halmd::mdsim::gpu::particle_kernel;

namespace halmd
{
namespace observables { namespace gpu
{
namespace profiles_kernel
{

/** number of particles */
__constant__ uint nbox_;

/**
 * compute bin indices from particle positions
 *
 * @return g_bin contains a bin ID (or offset) for each particle,
 * g_index holds the corresponding particle offset in g_r
 *
 * @param g_r tagged particle positions
 *
 * @param g_bin return value.
 * @param g_index return value.
 *
 * g_r, g_bin, and g_index must point to arrays of nbox_ elements
 */
template <typename vector_type, typename index_type>
__global__ void compute_bins(
    float4 const* g_r, uint* g_bin, uint* g_index
  , index_type ngrid, vector_type spacing, vector_type box_origin
)
{
    enum { dimension = vector_type::static_size };

    if (GTID >= nbox_) { return; }

    // read particle position
    vector_type r;
    unsigned type;
    tie(r, type) = untagged<vector_type>(g_r[GTID]);

    index_type index = element_mod(
        static_cast<index_type>(element_div(r - box_origin, spacing))
      , ngrid
    );

    // convert bin index to global bin ID
    // e.g., ID = x + ngrid_[0] * (y + ngrid_[1] * z)
    unsigned id = index[dimension - 1];
    for (int i = dimension - 2; i >= 0; i--) {
        id *= ngrid[i];
        id += index[i];
    }
    g_bin[GTID] = id;

    // store particle index
    if (g_index) {
        g_index[GTID] = GTID;
    }
}

/**
 * find boundaries in sorted list of bin IDs
 *
 * @return g_boundaries contains begin and end index for each bin,
 * it should be initialised (with a negative value) before
 * since empty bins are left untouched.
 *
 * @param g_bin sorted list of bin IDs.
 * Must point to an array of nbox_ elements.
 *
 * @param g_boundaries return value.
 * The size of g_boundaries must comply with the numbers of bins,
 * i.e., it must not be exceeded by the largest bin ID stored in g_bin.
 *
 * each thread reads the bin IDs from two adjacent indices
 * and stores these indices if the IDs differ
 */
__global__ void find_boundaries(uint const* g_bin, int2* g_boundaries)
{
    if (GTID >= nbox_) { return; }

    const unsigned int j = g_bin[GTID];
    // first index marks the start of a bin
    if (GTID == 0) {
        g_boundaries[j].x = GTID;
        return;
    }

    const unsigned int k = g_bin[GTID - 1];
    if (k < j) {
        // index marks the start of a bin
        g_boundaries[j].x = GTID;
        // previous index marks the end of another bin,
        // index points behind the range analogous to end iterator
        g_boundaries[k].y = GTID;
    }

    // last index marks the end of a bin
    if (GTID == nbox_ - 1) {
        g_boundaries[j].y = nbox_;  //< denotes the global end of the particle range
    }
}

/**
  * collect stress tensors from particles and compute
  * diagonal elements of total stress tensor for each bin
  *
  * @param g_boundaries offset ranges in g_index for each bin
  * @param g_index particle index for each offset
  * @param g_stress_pot_first potential part of stress tensor for each particle
  *        (from force module, first part only)
  * @param g_v particle velocities (required for kinetic part of stress tensor)
  * @param g_stress_diag return value
  *
  * g_boundaries and g_stress_diag must point at arrays of size prod(ngrid_),
  * g_index, g_stress_pot_first, and g_v must point at arrays of nbox_ values
  */
template <
    typename vector_type
  , typename gpu_stress_tensor_first_type
  , typename coalesced_vector_type
>
__global__ void collect_stress_tensor(
    int2 const* g_boundaries
  , uint const* g_index
  , gpu_stress_tensor_first_type const* g_stress_pot_first
  , float4 const* g_v
  , coalesced_vector_type* g_stress_diag
  , uint nbins
)
{
    enum { dimension = vector_type::static_size };

    if (GTID >= nbins) { return; }

    vector_type stress_diag(0);

    // loop over all particles in bin
    int2 bin = g_boundaries[GTID];
    for (int offset = bin.x; offset < bin.y; ++offset) {
        uint i = g_index[offset];            //< determine particle index

        // add diagonal elements of potential part
        stress_diag += static_cast<vector_type>(g_stress_pot_first[i]);

        // construct kinetic part of stress tensor
        vector_type v;
        uint tag;
        // drop high-precision part in case of dsfloat type
        tie(v, tag) = untagged<vector_type>(g_v[i]);
        // make kinetic part and add diagonal elements only
        stress_diag += static_cast<vector_type>(get<0>(split(mdsim::make_stress_tensor(v))));
    }

    // write result to global memory
    g_stress_diag[GTID] = stress_diag;
}


} // namespace profiles_kernel

template <int dimension>
profiles_wrapper<dimension> const profiles_wrapper<dimension>::kernel = {
    profiles_kernel::nbox_
  , profiles_kernel::compute_bins<vector_type, index_type>
  , profiles_kernel::find_boundaries
  , profiles_kernel::collect_stress_tensor<vector_type>
};

template class profiles_wrapper<3>;
template class profiles_wrapper<2>;

}} //namespace mdsim::gpu

} //namespace halmd
