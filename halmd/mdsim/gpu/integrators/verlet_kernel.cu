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
template <int dimension, typename float_type, typename gpu_vector_type>
__global__ void integrate(
    float4* g_position
  , gpu_vector_type* g_image
  , float4* g_velocity
  , gpu_vector_type const* g_force
  , float timestep
  , fixed_vector<float, dimension> box_length
)
{
    typedef fixed_vector<float_type, dimension> vector_type;
    typedef fixed_vector<float, dimension> float_vector_type;

    // kernel execution parameters
    unsigned int const thread = GTID;
    unsigned int const nthread = GTDIM;

    // read position, species, velocity, mass, image, force from global memory
    vector_type r, v;
    unsigned int species;
    float mass;
#ifdef USE_VERLET_DSFUN
    tie(r, species) <<= tie(g_position[thread], g_position[thread + nthread]);
    tie(v, mass) <<= tie(g_velocity[thread], g_velocity[thread + nthread]);
#else
    tie(r, species) <<= g_position[thread];
    tie(v, mass) <<= g_velocity[thread];
#endif
    float_vector_type f = g_force[thread];

    // advance position by full step, velocity by half step
    v += f * (timestep / 2) / mass;
    r += v * timestep;
    float_vector_type image = box_kernel::reduce_periodic(r, box_length);

    // store position, species, velocity, mass, image in global memory
#ifdef USE_VERLET_DSFUN
    tie(g_position[thread], g_position[thread + nthread]) <<= tie(r, species);
    tie(g_velocity[thread], g_velocity[thread + nthread]) <<= tie(v, mass);
#else
    g_position[thread] <<= tie(r, species);
    g_velocity[thread] <<= tie(v, mass);
#endif
    if (!(image == float_vector_type(0))) {
        g_image[thread] = image + static_cast<float_vector_type>(g_image[thread]);
    }
}

/**
 * Second leapfrog half-step of velocity-Verlet algorithm
 */
template <int dimension, typename float_type, typename gpu_vector_type>
__global__ void finalize(
    float4* g_velocity
  , gpu_vector_type const* g_force
  , float timestep
)
{
    // kernel execution parameters
    unsigned int const thread = GTID;
    unsigned int const nthread = GTDIM;

    // read velocity, mass, force from global memory
    fixed_vector<float_type, dimension> v;
    float mass;
#ifdef USE_VERLET_DSFUN
    tie(v, mass) <<= tie(g_velocity[thread], g_velocity[thread + nthread]);
#else
    tie(v, mass) <<= g_velocity[thread];
#endif
    fixed_vector<float, dimension> f = g_force[thread];

    // advance velocity by half step
    v += f * (timestep / 2) / mass;

    // store velocity, mass in global memory
#ifdef USE_VERLET_DSFUN
    tie(g_velocity[thread], g_velocity[thread + nthread]) <<= tie(v, mass);
#else
    g_velocity[thread] <<= tie(v, mass);
#endif
}

} // namespace verlet_kernel

template <int dimension>
verlet_wrapper<dimension> const verlet_wrapper<dimension>::wrapper = {
#ifdef USE_VERLET_DSFUN
    verlet_kernel::integrate<dimension, dsfloat>
  , verlet_kernel::finalize<dimension, dsfloat>
#else
    verlet_kernel::integrate<dimension, float>
  , verlet_kernel::finalize<dimension, float>
#endif
};

template class verlet_wrapper<3>;
template class verlet_wrapper<2>;

} // namespace mdsim
} // namespace gpu
} // namespace integrators
} // namespace halmd
