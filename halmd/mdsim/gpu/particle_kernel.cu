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
#include <halmd/utility/gpu/texfetch.cuh>
#include <halmd/utility/tuple.hpp>
#include "particle_kernel.hpp"

namespace halmd {
namespace mdsim {
namespace gpu {
namespace particle_kernel {

/** number of particles in simulation box */
static __constant__ unsigned int nbox_;
/** number of particle types */
static __constant__ unsigned int ntype_;

/**
 * rearrange particles by a given permutation
 */
template <int dimension, typename float_type, typename ptr_type, typename aligned_vector_type>
__global__ void rearrange(
    cudaTextureObject_t t_r
  , cudaTextureObject_t t_image
  , cudaTextureObject_t t_v
  , cudaTextureObject_t t_id
  , unsigned int const* g_index
  , ptr_type g_r
  , aligned_vector_type* g_image
  , ptr_type g_v
  , unsigned int* g_id
  , unsigned int npart
)
{
    typedef fixed_vector<float_type, dimension> vector_type;
    int const i = (GTID < npart) ? g_index[GTID] : GTID;

    // copy position and velocity as float4 values, and image vector
    g_r[GTID] = texFetch<float4, float_type>::fetch(t_r, i);
    g_v[GTID] = texFetch<float4, float_type>::fetch(t_v, i);

    // select correct image texture depending on the space dimension
    g_image[GTID] = tex1Dfetch<aligned_vector_type>(t_image, i);

    // copy particle IDs
    g_id[GTID] = tex1Dfetch<unsigned int>(t_id, i);
}

} // namespace particle_kernel

template <int dimension, typename float_type>
particle_wrapper<dimension, float_type> particle_wrapper<dimension, float_type>::kernel = {
    particle_kernel::nbox_
  , particle_kernel::ntype_
  , particle_kernel::rearrange<dimension, float_type, ptr_type>
};

template class particle_wrapper<3, float>;
template class particle_wrapper<2, float>;
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
template class particle_wrapper<3, dsfloat>;
template class particle_wrapper<2, dsfloat>;
#endif

template<typename T>
static __global__ void particle_initialize_kernel (
  T *g_v
, T value
, T ghost_value
, unsigned int nparticle
)
{
    g_v[GTID] = (GTID < nparticle) ? value : ghost_value;
}

template<typename T>
particle_initialize_wrapper<T> particle_initialize_wrapper<T>::kernel = {
  particle_initialize_kernel<T>
};

#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
template<typename ptr_type, typename type>
static __global__ void dsfloat_particle_initialize_kernel (
  ptr_type g_v
  , type value
  , type ghost_value
  , unsigned int nparticle
)
{
    g_v[GTID] = make_tuple((GTID < nparticle) ? value : ghost_value, type());
}

template<size_t N>
dsfloat_particle_initialize_wrapper<N> dsfloat_particle_initialize_wrapper<N>::kernel = {
  dsfloat_particle_initialize_kernel<ptr_type, type>
};
#endif // USE_GPU_DOUBLE_SINGLE_PRECISION

template class particle_initialize_wrapper<unsigned int>;
template class particle_initialize_wrapper<uint2>;
template class particle_initialize_wrapper<uint4>;
template class particle_initialize_wrapper<int>;
template class particle_initialize_wrapper<int2>;
template class particle_initialize_wrapper<int4>;
template class particle_initialize_wrapper<float>;
template class particle_initialize_wrapper<float2>;
template class particle_initialize_wrapper<float4>;
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
template class dsfloat_particle_initialize_wrapper<1>;
template class dsfloat_particle_initialize_wrapper<2>;
template class dsfloat_particle_initialize_wrapper<4>;
#endif

} // namespace gpu
} // namespace mdsim
} // namespace halmd
