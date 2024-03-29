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

namespace halmd {
namespace mdsim {
namespace gpu {
namespace particle_group_kernel {

/**
 * copy a subset of a particle instance (particle group) to another particle
 * instance
 */
template <typename float_type, typename vector_type, typename ptr_type,
    typename aligned_vector_type>
__global__ void particle_group_to_particle(
    cudaTextureObject_t r
  , cudaTextureObject_t image
  , cudaTextureObject_t v
  , unsigned int const* g_index
  , ptr_type g_r
  , aligned_vector_type* g_image
  , ptr_type g_v
  , unsigned int npart
)
{
    enum { dimension = vector_type::static_size };

    if (GTID < npart) {
        int const i = g_index[GTID];

        // copy position and velocity as float4 values, and image vector
        g_r[GTID] = texFetch<float4, float_type>::fetch(r, i);
        g_v[GTID] = texFetch<float4, float_type>::fetch(v, i);

        // copy image vector with its type depending on the space dimension
        //g_image[GTID] = tex1Dfetch<typename particle_group_wrapper<dimension, float>::aligned_vector_type>(image, i);
        g_image[GTID] = tex1Dfetch<aligned_vector_type>(image, i);
    }
}

} // namespace particle_group_kernel

template <int dimension, typename float_type>
particle_group_wrapper<dimension, float_type>
particle_group_wrapper<dimension, float_type>::kernel = {
    particle_group_kernel::particle_group_to_particle<float_type, fixed_vector<float_type, dimension>, ptr_type>
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
