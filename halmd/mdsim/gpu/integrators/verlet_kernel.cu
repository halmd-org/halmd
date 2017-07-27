/*
 * Copyright © 2014      Felix Höfling
 * Copyright © 2008-2012 Peter Colberg
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

#include <halmd/mdsim/gpu/box_kernel.cuh>
#include <halmd/mdsim/gpu/integrators/verlet_kernel.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>
#include <halmd/utility/gpu/thread.cuh>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace integrators {
namespace verlet_kernel {

/**
 * First leapfrog half-step of velocity-Verlet algorithm
 */
template <int dimension, typename float_type, typename ptr_type, typename gpu_vector_type>
__global__ void integrate(
    ptr_type g_position
  , gpu_vector_type* g_image
  , ptr_type g_velocity
  , gpu_vector_type const* g_force
  , float timestep
  , fixed_vector<float, dimension> box_length
)
{
    typedef fixed_vector<float_type, dimension> vector_type;
    typedef fixed_vector<float, dimension> float_vector_type;

    // kernel execution parameters
    unsigned int const thread = GTID;

    // read position, species, velocity, mass, image, force from global memory
    vector_type r, v;
    unsigned int species;
    float mass;
    tie(r, species) <<= g_position[thread];
    tie(v, mass) <<= g_velocity[thread];
    float_vector_type f = g_force[thread];

    // advance position by full step, velocity by half step
    v += f * (timestep / 2) / mass;
    r += v * timestep;
    float_vector_type image = box_kernel::reduce_periodic(r, box_length);

    // store position, species, velocity, mass, image in global memory
    g_position[thread] <<= tie(r, species);
    g_velocity[thread] <<= tie(v, mass);
    if (!(image == float_vector_type(0))) {
        g_image[thread] = image + static_cast<float_vector_type>(g_image[thread]);
    }
}

/**
 * Second leapfrog half-step of velocity-Verlet algorithm
 */
template <int dimension, typename float_type, typename ptr_type, typename gpu_vector_type>
__global__ void finalize(
    ptr_type g_velocity
  , gpu_vector_type const* g_force
  , float timestep
)
{
    // kernel execution parameters
    unsigned int const thread = GTID;

    // read velocity, mass, force from global memory
    fixed_vector<float_type, dimension> v;
    float mass;
    tie(v, mass) <<= g_velocity[thread];
    fixed_vector<float, dimension> f = g_force[thread];

    // advance velocity by half step
    v += f * (timestep / 2) / mass;

    // store velocity, mass in global memory
    g_velocity[thread] <<= tie(v, mass);
}

} // namespace verlet_kernel

template <int dimension, typename float_type>
verlet_wrapper<dimension, float_type> const verlet_wrapper<dimension, float_type>::wrapper = {
    verlet_kernel::integrate<dimension, float_type, ptr_type>
  , verlet_kernel::finalize<dimension, float_type, ptr_type>
};

#ifdef USE_GPU_SINGLE_PRECISION
template class verlet_wrapper<3, float>;
template class verlet_wrapper<2, float>;
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
template class verlet_wrapper<3, dsfloat>;
template class verlet_wrapper<2, dsfloat>;
#endif

} // namespace mdsim
} // namespace gpu
} // namespace integrators
} // namespace halmd
