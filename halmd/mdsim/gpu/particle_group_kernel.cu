/*
 * Copyright © 2014 Nicolas Höft
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

#include <halmd/mdsim/gpu/particle_group_kernel.hpp>
#include <halmd/utility/gpu/thread.cuh>
#include <halmd/utility/gpu/texfetch.cuh>
#include <halmd/utility/tuple.hpp>

/** positions, types */
static texture<float4> r_;
/** orientations */
static texture<float4> u_;
/** velocities, masses */
static texture<float4> v_;

namespace halmd {
namespace mdsim {
namespace gpu {
namespace particle_group_kernel {

/** positions, types */
static texture<float4> r_;
/** orientations */
static texture<float4> u_;
/** velocities, masses */
static texture<float4> v_;

/** minimum image vectors */
template<int dimension>
struct image
{
    // instantiate a separate texture for each aligned vector type
    typedef texture<typename particle_group_wrapper<dimension, float>::aligned_vector_type> type;
    static type tex_;
};
// instantiate static members
template<int dimension> image<dimension>::type image<dimension>::tex_;

/**
 * copy a subset of a particle instance (particle group) to another particle instance
 */
template <typename float_type, typename vector_type, typename ptr_type, typename aligned_vector_type>
__global__ void particle_group_to_particle(
    unsigned int const* g_index
  , ptr_type g_r
  , aligned_vector_type* g_image
  , ptr_type g_u
  , ptr_type g_v
  , unsigned int npart
)
{
    enum { dimension = vector_type::static_size };

    if (GTID < npart) {
        int const i = g_index[GTID];

        // copy position and velocity as float4 values, and image vector
        g_r[GTID] = texFetch<float_type>::fetch(r_, i);
        g_u[GTID] = texFetch<float_type>::fetch(u_, i);
        g_v[GTID] = texFetch<float_type>::fetch(v_, i);

        // copy image vector with its type depending on the space dimension
        g_image[GTID] = tex1Dfetch(image<dimension>::tex_, i);
    }
}

} // namespace particle_group_kernel

template <int dimension, typename float_type>
particle_group_wrapper<dimension, float_type> const
particle_group_wrapper<dimension, float_type>::kernel = {
    particle_group_kernel::r_
  , particle_group_kernel::image<dimension>::tex_
  , particle_group_kernel::u_
  , particle_group_kernel::v_
  , particle_group_kernel::particle_group_to_particle<float_type, fixed_vector<float_type, dimension>, ptr_type>
};

template class particle_group_wrapper<3, float>;
template class particle_group_wrapper<2, float>;
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
template class particle_group_wrapper<3, dsfloat>;
template class particle_group_wrapper<2, dsfloat>;
#endif

} // namespace gpu
} // namespace mdsim
} // namespace halmd
