/*
 * Copyright © 2014 Nicolas Höft
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

#include <halmd/mdsim/gpu/particle_group_kernel.hpp>
#include <halmd/utility/gpu/thread.cuh>
#include <halmd/utility/tuple.hpp>

/** positions, types */
static texture<float4> r_;
/** velocities, masses */
static texture<float4> v_;

namespace halmd {
namespace mdsim {
namespace gpu {
namespace particle_group_kernel {

/** positions, types */
static texture<float4> r_;
/** velocities, masses */
static texture<float4> v_;

/** minimum image vectors */
template<int dimension>
struct image
{
    // instantiate a separate texture for each aligned vector type
    typedef texture<typename particle_group_wrapper<dimension>::aligned_vector_type> type;
    static type tex_;
};
// instantiate static members
template<int dimension> image<dimension>::type image<dimension>::tex_;

/**
 * copy a subset of a particle instance (particle group) to another particle instance
 */
template <typename vector_type, typename aligned_vector_type>
__global__ void particle_group_to_particle(
    unsigned int const* g_index
  , float4* g_r
  , aligned_vector_type* g_image
  , float4* g_v
  , unsigned int npart
)
{
    enum { dimension = vector_type::static_size };

    if (GTID < npart) {
        int const i = g_index[GTID];

        // copy position and velocity as float4 values, and image vector
        g_r[GTID] = tex1Dfetch(r_, i);
        g_v[GTID] = tex1Dfetch(v_, i);
#ifdef USE_VERLET_DSFUN
        g_r[GTID + GTDIM] = tex1Dfetch(r_, i + GTDIM);
        g_v[GTID + GTDIM] = tex1Dfetch(v_, i + GTDIM);
#endif

        // copy image vector with its type depending on the space dimension
        g_image[GTID] = tex1Dfetch(image<dimension>::tex_, i);
    }
}

} // namespace particle_group_kernel

template <int dimension>
particle_group_wrapper<dimension> const
particle_group_wrapper<dimension>::kernel = {
    particle_group_kernel::r_
  , particle_group_kernel::image<dimension>::tex_
  , particle_group_kernel::v_
#ifdef USE_VERLET_DSFUN
  , particle_group_kernel::particle_group_to_particle<fixed_vector<dsfloat, dimension> >
#else
  , particle_group_kernel::particle_group_to_particle<fixed_vector<float, dimension> >
#endif
};

template class particle_group_wrapper<3>;
template class particle_group_wrapper<2>;

} // namespace gpu
} // namespace mdsim
} // namespace halmd
