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
#include <halmd/mdsim/gpu/particle_kernel.cuh>
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
/** ids */
static texture<unsigned int> id_;

/** minimum image vectors */
template<int dimension>
struct image
{
    // instantiate a separate texture for each aligned vector type
    typedef texture<typename particle_wrapper<float, dimension>::aligned_vector_type> type;
    static type tex_;
};
// instantiate static members
template<int dimension> image<dimension>::type image<dimension>::tex_;

/**
 * initialize particle positions and species, velocity and mass
 */
 template<typename float_type>
__global__ void initialize(
    typename type_traits<4, float_type>::gpu::ptr_type g_r
  , typename type_traits<4, float_type>::gpu::ptr_type g_v
  , unsigned int size
)
{
    unsigned int type = (GTID < size) ? 0 : placeholder;
    fixed_vector<float_type, 3> r (0.0f);
    g_r[GTID] <<= tie(r, type);
    g_v[GTID] <<= make_tuple(r, 1.0f);
}

/**
 * rearrange particles by a given permutation
 */
template <typename float_type, int dimension, typename aligned_vector_type>
__global__ void rearrange(
    unsigned int const* g_index
  , typename type_traits<4, float_type>::gpu::ptr_type g_r
  , aligned_vector_type* g_image
  , typename type_traits<4, float_type>::gpu::ptr_type g_v
  , unsigned int* g_id
  , unsigned int npart
)
{
    typedef fixed_vector<float_type, dimension> vector_type;
    /*if (GTID < npart)*/ {
        int const i = (GTID < npart) ? g_index[GTID] : GTID;

        // copy position and velocity as float4 values, and image vector
        g_r[GTID] = texFetch<float_type>::fetch(r_, i);
        g_v[GTID] = texFetch<float_type>::fetch(v_, i);

        // select correct image texture depending on the space dimension
        g_image[GTID] = tex1Dfetch(image<dimension>::tex_, i);

        // copy particle ids
        g_id[GTID] = tex1Dfetch(id_, i);
    }
}

} // namespace particle_kernel

template <typename float_type, int dimension>
particle_wrapper<float_type, dimension> const particle_wrapper<float_type, dimension>::kernel = {
    particle_kernel::nbox_
  , particle_kernel::ntype_
  , particle_kernel::ntypes_
  , particle_kernel::r_
  , particle_kernel::image<dimension>::tex_
  , particle_kernel::v_
  , particle_kernel::id_
  , particle_kernel::initialize<float_type>
  , particle_kernel::rearrange<float_type, dimension>
};

template class particle_wrapper<float, 3>;
template class particle_wrapper<float, 2>;
template class particle_wrapper<dsfloat, 3>;
template class particle_wrapper<dsfloat, 2>;

} // namespace gpu
} // namespace mdsim
} // namespace halmd
