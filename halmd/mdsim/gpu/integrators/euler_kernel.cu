/*
 * Copyright © 2011-2014 Felix Höfling
 * Copyright © 2011-2012 Michael Kopp
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
#include <halmd/mdsim/gpu/integrators/euler_kernel.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>
#include <halmd/utility/gpu/thread.cuh>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace integrators {
namespace euler_kernel {

/**
 * Euler integration: @f$ r(t + \Delta t) = r(t) + v(t) \Delta t @f$
 *
 * @param g_position    positions
 * @param g_image       number of times the particle exceeded the box margin
 * @param g_velocity    velocities
 * @param g_force       forces
 * @param timestep      integration timestep
 * @param box_length    edge lengths of cuboid box
 */
template <int dimension, typename float_type, typename gpu_vector_type>
__global__ void integrate(
    float4* g_position
  , gpu_vector_type* g_image
  , float4 const* g_velocity
  , float timestep
  , fixed_vector<float, dimension> box_length
)
{
    typedef fixed_vector<float_type, dimension> vector_type;
    typedef fixed_vector<float, dimension> float_vector_type;

    // kernel execution parameters
    unsigned int const thread = GTID;
    unsigned int const nthread = GTDIM;

    // read position, species, velocity, mass, image from global memory
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

    // Euler integration
    r += v * timestep;
    // enforce periodic boundary conditions
    float_vector_type image = box_kernel::reduce_periodic(r, box_length);

    // store position, species, image in global memory
#ifdef USE_VERLET_DSFUN
    tie(g_position[thread], g_position[thread + nthread]) <<= tie(r, species);
#else
    g_position[thread] <<= tie(r, species);
#endif
    if (!(image == float_vector_type(0))) {
        g_image[thread] = image + static_cast<float_vector_type>(g_image[thread]);
    }
}

} // namespace euler_kernel

template <int dimension>
euler_wrapper<dimension> const euler_wrapper<dimension>::kernel = {
#ifdef USE_VERLET_DSFUN
    euler_kernel::integrate<dimension, dsfloat>
#else
    euler_kernel::integrate<dimension, float>
#endif
};

// explicit instantiation
template class euler_wrapper<3>;
template class euler_wrapper<2>;

} // namespace integrators
} // namespace gpu
} // namespace mdsim
} // namespace halmd
