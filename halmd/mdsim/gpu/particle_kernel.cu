/*
 * Copyright © 2010-2011  Peter Colberg and Felix Höfling
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
#include <halmd/utility/gpu/variant.cuh>

using namespace halmd::algorithm::gpu;
using namespace halmd::utility::gpu;

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
/** positions, types */
static texture<float4> r_;
/** minimum image vectors */
static texture<variant<map<pair<int_<3>, float4>, pair<int_<2>, float2> > > > image_;
/** velocities, tags */
static texture<float4> v_;

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

/**
 * generate ascending index sequence
 */
__global__ void gen_index(unsigned int* g_index)
{
    g_index[GTID] = (GTID < nbox_) ? GTID : 0;
}

/**
 * rearrange particles by a given permutation
 */
template <typename vector_type, typename aligned_vector_type>
__global__ void rearrange(
    unsigned int const* g_index
  , float4* g_r
  , aligned_vector_type* g_image
  , float4* g_v
  , unsigned int* g_tag
)
{
    enum { dimension = vector_type::static_size };

    int const i = g_index[GTID];

    // copy position including type, and image vector
    g_r[GTID] = tex1Dfetch(r_, i);
#ifdef USE_VERLET_DSFUN
    g_r[GTID + GTDIM] = tex1Dfetch(r_, i + GTDIM);
#endif
    g_image[GTID] = tex1Dfetch(get<dimension>(image_), i);

    // copy velocity, but split off tag and store separately
    {
        vector_type v;
        unsigned int tag;
#ifdef USE_VERLET_DSFUN
        tie(v, tag) = untagged<vector_type>(tex1Dfetch(v_, i), tex1Dfetch(v_, i + GTDIM));
        tie(g_v[GTID], g_v[GTID + GTDIM]) = tagged(v, tag);
#else
        tie(v, tag) = untagged<vector_type>(tex1Dfetch(v_, i));
        g_v[GTID] = tagged(v, tag);
#endif
        g_tag[GTID] = tag;
    }
}

} // namespace particle_kernel

template <int dimension>
particle_wrapper<dimension> const particle_wrapper<dimension>::kernel = {
    particle_kernel::nbox_
  , particle_kernel::ntype_
  , particle_kernel::ntypes_
  , particle_kernel::r_
  , get<dimension>(particle_kernel::image_)
  , particle_kernel::v_
  , particle_kernel::tag<fixed_vector<float, dimension> >
  , particle_kernel::gen_index
#ifdef USE_VERLET_DSFUN
  , particle_kernel::rearrange<fixed_vector<dsfloat, dimension> >
#else
  , particle_kernel::rearrange<fixed_vector<float, dimension> >
#endif
};

template class particle_wrapper<3>;
template class particle_wrapper<2>;

} // namespace mdsim
} // namespace gpu
} // namespace halmd
