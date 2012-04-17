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

#include <halmd/mdsim/gpu/box_kernel.cuh>
#include <halmd/mdsim/gpu/particle_kernel.cuh>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/observables/gpu/phase_space_kernel.hpp>
#include <halmd/utility/gpu/thread.cuh>

using namespace halmd::mdsim::gpu; //< namespace box_kernel

namespace halmd {
namespace observables {
namespace gpu {
namespace phase_space_kernel {

/** positions, types */
texture<float4> r_;
/** minimum image vectors */
texture<void> image_;
/** velocities, tags */
texture<float4> v_;

/**
 * sample phase space for all particle of a single species
 */
template <typename vector_type, typename T>
__global__ void sample(
    unsigned int const* g_reverse_tag
  , T* g_r
  , T* g_v
  , vector_type box_length
  , unsigned int npart
)
{
    using mdsim::gpu::particle_kernel::untagged;
    using mdsim::gpu::particle_kernel::tagged;

    enum { dimension = vector_type::static_size };
    typedef typename phase_space_wrapper<dimension>::coalesced_vector_type coalesced_vector_type;

    if (GTID < npart) {
        // permutation index
        uint const rtag = g_reverse_tag[GTID];
        // fetch particle from texture caches
        unsigned int tag, type;
        vector_type r, v;
        tie(r, type) = untagged<vector_type>(tex1Dfetch(r_, rtag));
        tie(v, tag) = untagged<vector_type>(tex1Dfetch(v_, rtag));
        // extend particle positions in periodic box
        vector_type image = tex1Dfetch(reinterpret_cast<texture<coalesced_vector_type>&>(image_), rtag);
        box_kernel::extend_periodic(r, image, box_length);
        // store particle in global memory
        g_r[GTID] = tagged(r, type);
        g_v[GTID] = tagged(v, type);
    }
}

/**
 * shift particle positions to range (-L/2, L/2)
 */
template <typename vector_type, typename coalesced_vector_type>
__global__ void reduce_periodic(
    unsigned int const* g_reverse_tag
  , float4* g_r
  , coalesced_vector_type* g_image
  , vector_type box_length
  , unsigned int npart
)
{
    using mdsim::gpu::particle_kernel::untagged;
    using mdsim::gpu::particle_kernel::tagged;

    enum { dimension = vector_type::static_size };

    if (GTID < npart) {
        unsigned int rtag = g_reverse_tag[GTID];
        vector_type r;
        unsigned int type;
        tie(r, type) = untagged<vector_type>(tex1Dfetch(r_, rtag));

        vector_type image = box_kernel::reduce_periodic(r, box_length);

        g_image[rtag] = image;
        g_r[rtag] = tagged(r, type);
    }
}

} // namespace phase_space_kernel

template <int dimension>
phase_space_wrapper<dimension> const phase_space_wrapper<dimension>::kernel = {
    phase_space_kernel::r_
  , phase_space_kernel::image_
  , phase_space_kernel::v_
  , phase_space_kernel::sample<fixed_vector<float, dimension> >
  , phase_space_kernel::reduce_periodic<fixed_vector<float, dimension> >
};

template class phase_space_wrapper<3>;
template class phase_space_wrapper<2>;

} // namespace observables
} // namespace gpu
} // namespace halmd
