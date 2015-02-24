/*
 * Copyright © 2008-2010  Peter Colberg
 * Copyright © 2015       Nicolas Höft
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
/** velocities, tags */
texture<float4> v_;

/** minimum image vectors */
template<int dimension>
struct image
{
    // instantiate a separate texture for each aligned vector type
    typedef texture<typename phase_space_wrapper<dimension>::coalesced_vector_type> type;
    static type tex_;
};
// instantiate static members
template<int dimension> image<dimension>::type image<dimension>::tex_;

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
    enum { dimension = vector_type::static_size };
    typedef typename phase_space_wrapper<dimension>::coalesced_vector_type coalesced_vector_type;

    if (GTID < npart) {
        // permutation index
        uint const rtag = g_reverse_tag[GTID];
        // fetch particle from texture caches
        unsigned int tag, type;
        vector_type r, v;
        tie(r, type) <<= tex1Dfetch(r_, rtag);
        tie(v, tag) <<= tex1Dfetch(v_, rtag);
        // extend particle positions in periodic box
        vector_type img = tex1Dfetch(image<dimension>::tex_, rtag);
        box_kernel::extend_periodic(r, img, box_length);
        // store particle in global memory
        g_r[GTID] <<= tie(r, type);
        g_v[GTID] <<= tie(v, type);
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
    enum { dimension = vector_type::static_size };

    if (GTID < npart) {
        unsigned int rtag = g_reverse_tag[GTID];
        vector_type r;
        unsigned int type;
        tie(r, type) <<= tex1Dfetch(r_, rtag);

        vector_type image = box_kernel::reduce_periodic(r, box_length);

        g_image[rtag] = image;
        g_r[rtag] <<= tie(r, type);
    }
}

/**
 * copy particle group from GPU to host memory
 */
__global__ void copy_particle_group(
    unsigned int const* g_group
  , unsigned int* h_group
  , unsigned int size
)
{
    if (GTID < size) {
        h_group[GTID] = g_group[GTID];
    }
}

} // namespace phase_space_kernel

template <int dimension>
phase_space_wrapper<dimension> const phase_space_wrapper<dimension>::kernel = {
    phase_space_kernel::r_
  , phase_space_kernel::image<dimension>::tex_
  , phase_space_kernel::v_
  , phase_space_kernel::sample<fixed_vector<float, dimension> >
  , phase_space_kernel::reduce_periodic<fixed_vector<float, dimension> >
  , phase_space_kernel::copy_particle_group
};

template class phase_space_wrapper<3>;
template class phase_space_wrapper<2>;

} // namespace gpu
} // namespace observables
} // namespace halmd
