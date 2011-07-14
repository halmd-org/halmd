/*
 * Copyright © 2010  Peter Colberg
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

#include <halmd/algorithm/gpu/tuple.cuh>
#include <halmd/mdsim/gpu/particle_kernel.cuh>
#include <halmd/mdsim/gpu/particle_kernel.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/utility/gpu/thread.cuh>

using namespace halmd::algorithm::gpu;

namespace halmd {
namespace mdsim {
namespace gpu {
namespace particle_kernel {

/** number of particles in simulation box */
static __constant__ unsigned int nbox_;
/** number of particle types */
static __constant__ unsigned int ntype_;
/** number of particles per type */
static texture<unsigned int> ntypes_;

/**
 * set particle tags and types
 */
template <
    typename vector_type
  , typename coalesced_vector_type
>
__global__ void tag(coalesced_vector_type* g_r, coalesced_vector_type* g_v)
{
    vector_type r, v;
    unsigned int type, tag;
    tie(r, type) = untagged<vector_type>(g_r[GTID]);
    tie(v, tag) = untagged<vector_type>(g_v[GTID]);

    // set particle identifier unique within each particle type,
    // use a 0-based continuous numbering
    tag = GTID;

    // set particle type and adjust tag
    for (type = 0; type < ntype_; ++type) {
        unsigned int n = tex1Dfetch(ntypes_, type);
        if (tag < n) {
            break;
        }
        tag -= n;
    }

    g_r[GTID] = tagged(r, type);
    g_v[GTID] = tagged(v, tag);
}

} // namespace particle_kernel

template <int dimension>
particle_wrapper<dimension> const particle_wrapper<dimension>::kernel = {
    particle_kernel::nbox_
  , particle_kernel::ntype_
  , particle_kernel::ntypes_
  , particle_kernel::tag<fixed_vector<float, dimension> >
};

template class particle_wrapper<3>;
template class particle_wrapper<2>;

} // namespace mdsim
} // namespace gpu
} // namespace halmd
