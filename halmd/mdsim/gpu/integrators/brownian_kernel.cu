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
#include <halmd/mdsim/gpu/integrators/brownian_kernel.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>
#include <halmd/random/gpu/normal_distribution.cuh>
#include <halmd/random/gpu/random_number_generator.cuh>
#include <halmd/utility/gpu/thread.cuh>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace integrators {
namespace brownian_kernel {

/**array of diffusion constants */
static texture<float> param_;

/**
 * Brownian integration: @f$ r(t + \Delta t) = r(t) + v(t) \Delta t @f$
 *
 * @param g_position    positions
 * @param g_image       number of times the particle exceeded the box margin
 * @param g_velocity    velocities
 * @param g_force       forces
 * @param timestep      integration timestep
 * @param D             diffusion constant 
 * @param rng           instance of random number generator 
 * @param box_length    edge lengths of cuboid box
 */
template <int dimension, typename float_type, typename rng_type, typename gpu_vector_type>
__global__ void integrate(
    float4* g_position
  , gpu_vector_type* g_image
  , float4 const* g_velocity
  , float timestep
  , rng_type rng
  , unsigned int nparticle
  , unsigned int nplace 
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
    float D;
    float const mean = 0.f;
    
    //read random number generator state from global device memory
    typename rng_type::state_type state = rng[thread];

    //initialize random displacement
    vector_type dr;
    float_type cache;
    bool cached = false;


    for (uint i = thread; i < nparticle; i += nthread) {
#ifdef USE_VERLET_DSFUN
        tie(r, species) <<= tie(g_position[i], g_position[i + nplace]);
        tie(v, mass) <<= tie(g_velocity[i], g_velocity[i + nplace]);
#else
        tie(r, species) <<= g_position[i];
        tie(v, mass) <<= g_velocity[i];
#endif

        //normal distribution parameters
        D = 1; // tex1Dfetch(param_, species);
       float const sigma = sqrtf(2 * D * timestep);
   
        //compute random displacement
        tie(dr[0], dr[1]) =  random::gpu::normal(rng, state, mean, sigma);
        //if (dimension % 2 ) {
            //if ((cached = !cached)) {
                tie(dr[2], cache) = random::gpu::normal(rng, state, mean, sigma);
            //}
            //else{
            //    dr[dimension - 1] = cache;
            //}
        //}
        // Brownian integration
        r += dr; 
        // enforce periodic boundary conditions
        float_vector_type image =  float_vector_type(0);// box_kernel::reduce_periodic(r, box_length);

        // store position, species, image in global memory
#ifdef USE_VERLET_DSFUN
        tie(g_position[i], g_position[i + nplace]) <<= tie(r, species);
#else
        g_position[i] <<= tie(r, species);
#endif
        if (!(image == float_vector_type(0))) {
            g_image[i] = image + static_cast<float_vector_type>(g_image[i]);
        }
    }
    // store random number generator state in global device memory
    rng[thread] = state;
}

} // namespace brownian_kernel

template <int dimension, typename rng_type>
cuda::texture<float> brownian_wrapper<dimension, rng_type>::param = brownian_kernel::param_;

template <int dimension, typename rng_type>
brownian_wrapper<dimension, rng_type> const brownian_wrapper<dimension, rng_type>::kernel = {
#ifdef USE_VERLET_DSFUN
    brownian_kernel::integrate<dimension, dsfloat, rng_type>
#else
    brownian_kernel::integrate<dimension, float, rng_type>
#endif
};

// explicit instantiation
template class brownian_wrapper<3, random::gpu::rand48_rng>;
template class brownian_wrapper<2, random::gpu::rand48_rng>;

} // namespace integrators
} // namespace gpu
} // namespace mdsim
} // namespace halmd
