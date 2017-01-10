/*
 * Copyright © 2010-2011 Felix Höfling
 * Copyright © 2015      Nicolas Höft
 * Copyright © 2010-2011 Peter Colberg
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

#include <halmd/mdsim/gpu/particle_kernel.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/utility/gpu/thread.cuh>
#include <halmd/utility/tuple.hpp>

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
/** velocities, masses */
static texture<float4> v_;
/** IDs */
static texture<unsigned int> id_;

/** minimum image vectors */
template<int dimension>
struct image
{
    // instantiate a separate texture for each aligned vector type
    typedef texture<typename particle_wrapper<dimension>::aligned_vector_type> type;
    static type tex_;
};
// instantiate static members
template<int dimension> image<dimension>::type image<dimension>::tex_;

/**
 * rearrange particles by a given permutation
 */
template <typename vector_type, typename aligned_vector_type>
__global__ void rearrange(
    unsigned int const* g_index
  , float4* g_r
  , aligned_vector_type* g_image
  , float4* g_v
  , unsigned int* g_id
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

        // select correct image texture depending on the space dimension
        g_image[GTID] = tex1Dfetch(image<dimension>::tex_, i);

        // copy particle IDs
        g_id[GTID] = tex1Dfetch(id_, i);
    }
}

} // namespace particle_kernel

template <int dimension>
particle_wrapper<dimension> const particle_wrapper<dimension>::kernel = {
    particle_kernel::nbox_
  , particle_kernel::ntype_
  , particle_kernel::ntypes_
  , particle_kernel::r_
  , particle_kernel::image<dimension>::tex_
  , particle_kernel::v_
  , particle_kernel::id_
#ifdef USE_VERLET_DSFUN
  , particle_kernel::rearrange<fixed_vector<dsfloat, dimension> >
#else
  , particle_kernel::rearrange<fixed_vector<float, dimension> >
#endif
};

template class particle_wrapper<3>;
template class particle_wrapper<2>;

} // namespace gpu
} // namespace mdsim
} // namespace halmd
