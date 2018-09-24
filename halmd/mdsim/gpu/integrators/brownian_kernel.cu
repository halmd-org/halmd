/*
 * Copyright © 2016      Manuel Dibak
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
#include <limits>
#include <cassert>
#include <cmath>
#include <cstdio>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace integrators {
namespace brownian_kernel {

/**array of diffusion constants */
static texture<float4> param_;

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

template <int dimension, typename float_type, typename vector_type>
__device__ void update_displacement(
    float_type const D_par
  , float_type const D_perp
  , float_type const prop_strength
  , vector_type & r
  , vector_type & u
  , vector_type & v
  , vector_type const& f
  , float timestep
  , float temperature
  , float_type eta1
  , float_type eta2
  , float_type eta3
)
{
    //random and systematic components of displacement
    vector_type dr_r, dr_s;

    dr_r[0] = eta1;
    dr_r[1] = eta2;
    if(dimension % 2) {
        dr_r[2] = eta3;
    }

    //now systematic part
    float_type f_u = inner_prod(f, u);
    dr_s = ( timestep * f_u * D_par / temperature + prop_strength * timestep)
        * u + ( timestep * D_perp / temperature) * (f - f_u * u);

    r += dr_r + dr_s;
}

template <typename float_type>
__device__ void update_orientation(
    float_type const D_rot
  , fixed_vector<float_type, 2> & u
  , fixed_vector<float_type, 2> const& tau
  , float timestep
  , float temp
  , float_type eta1
  , float_type eta2
  , float epsilon
  , float_type const sigma_rot
)
{

    fixed_vector<float_type, 2> e1;

    //construct trihedron along particle orientation for movement in x-y plane
    // e1 lies in x-y plane
    e1[0] = -u[1];
    e1[1] = u[0];

    fixed_vector<float_type, 1> omega;
    // first term is the random torque, second is the systematic
    // torque
    omega = eta1 + tau[0] * D_rot * timestep / temp;

    float  alpha = omega[0];
    u = cos( alpha ) * u + sin( alpha ) * e1;

    //ensure normalization
    if ( (float) norm_2(u) > 2e-38){
        u /= norm_2(u);
    }
}

template <typename float_type>
__device__ void update_orientation(
    float_type const D_rot
  , fixed_vector<float_type, 3> & u
  , fixed_vector<float_type, 3> const& tau
  , float timestep
  , float temp
  , float_type eta1
  , float_type eta2
  , float epsilon
  , float_type const sigma_rot
)
{
    fixed_vector<float_type, 3> e1, e2;

    //construct trihedron along particle orientation
    if ( (float) u[1] > epsilon || (float) u[2] > epsilon) {
        e1[0] = 0; e1[1] = u[2]; e1[2] = -u[1];
    }
    else {
        e1[0] = u[1]; e1[1] = -u[0]; e1[2] = 0;
    }
    e1 /= norm_2(e1);
    e2 = cross_prod(u, e1);
    e2 /= norm_2(e2);

    fixed_vector<float_type, 3> omega;

    // first two terms are the random angular velocity, the final is the
    // systematic torque
    omega = eta1 * e1 + eta2 * e2 + tau * D_rot * timestep / temp ;
    float  alpha = norm_2(omega);
    omega /= alpha;
    // Ω = eta1 * e1 + eta2 * e2
    // => Ω × u = (eta1 * e1 × u + eta2  * e2 × u) = eta2 * e1 - eta1 * e2
    u = (1 - cos(alpha)) * inner_prod(omega, u) * omega + cos(alpha) * u + sin(alpha) * cross_prod(omega, u);

    //ensure normalization
    if ( (float) norm_2(u) != 0){
        u /= norm_2(u);
    }
}

template <int dimension, typename float_type, typename rng_type, typename
    gpu_vector_type, typename gpu_pseudo_vector_type>
__global__ void integrate(
    float4* g_position
  , float4* g_orientation
  , gpu_vector_type* g_image
  , float4 const* g_velocity
  , gpu_vector_type const* g_force
  , gpu_pseudo_vector_type const* g_torque
  , float timestep
  , float temp
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

    //numerical limits for computation
    float epsilon = 1E-7;

    // read position, species, velocity, mass, image from global memory
    vector_type r, u, v;
    unsigned int species;
    float mass;
    float nothing;
    float const mean = 0;

    //read random number generator state from global device memory
    typename rng_type::state_type state = rng[thread];

    for (uint i = thread; i < nparticle; i += nthread) {

        // read in relevant variables
#ifdef USE_VERLET_DSFUN
        tie(r, species) <<= tie(g_position[i], g_position[i + nplace]);
        tie(u, nothing) <<= tie(g_orientation[i], g_orientation[i + nplace]);
        tie(v, mass)    <<= tie(g_velocity[i], g_velocity[i + nplace]);
#else
        tie(r, species) <<= g_position[i];
        tie(u, nothing) <<= g_orientation[i];
        tie(v, mass)    <<= g_velocity[i];
#endif

        //normal distribution parameters
        fixed_vector<float, 4> D    = tex1Dfetch(param_, species);
        float_type const D_perp     = D[0];
        float_type const D_par      = D[1];
        float_type const D_rot      = D[2];
        float_type const prop_str   = D[3];
        float_type const sigma_perp = sqrtf(2 * D_perp * timestep );
        float_type const sigma_par  = sqrtf(2 * D_par * timestep );
        float_type const sigma_rot  = sqrtf(2 * D_rot * timestep );

        vector_type f   = static_cast<float_vector_type>(g_force[i]);
        vector_type tau = static_cast<float_vector_type>(g_torque[i]);

        //draw random numbers
        float_type eta1, eta2, eta3, cache;
        bool cached = false;
        tie(eta1, eta2) =  random::gpu::normal(rng, state, mean, sigma_perp);
        if (dimension % 2 ) {
            if ((cached = !cached)) {
                tie(eta3, cache) = random::gpu::normal(rng, state, mean, sigma_par);
            }
            else{
                eta3 = cache;
            }
        }

        // Brownian integration
        update_displacement<dimension, float_type, vector_type>( D_par,
                D_perp, prop_str, r, u, v, f, timestep, temp, eta1, eta2,
                eta3);

        // enforce periodic boundary conditions
        float_vector_type image = box_kernel::reduce_periodic(r, box_length);


        tie(eta1, eta2) =  random::gpu::normal(rng, state, mean, sigma_rot);

        // update orientation last (Ito interpretation)
        update_orientation( D_rot, u, tau, timestep, temp, eta1, eta2,
                epsilon, sigma_rot);






        // store position, species, image in global memory
#ifdef USE_VERLET_DSFUN
        tie(g_position[i], g_position[i + nplace]) <<= tie(r, species);
        tie(g_orientation[i], g_orientation[i + nplace]) <<= tie(u, nothing);
#else
        g_position[i] <<= tie(r, species);
        g_orientation[i] <<= tie(u, nothing);
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
cuda::texture<float4> brownian_wrapper<dimension, rng_type>::param = brownian_kernel::param_;

template <int dimension, typename rng_type>
brownian_wrapper<dimension, rng_type> const brownian_wrapper<dimension, rng_type>::kernel = {
#ifdef USE_VERLET_DSFUN
    brownian_kernel::integrate<dimension, dsfloat, rng_type>
#else
    brownian_kernel::integrate<dimension, float, rng_type>
#endif
};

// explicit instantiation
template class brownian_wrapper<2, random::gpu::rand48_rng>;
template class brownian_wrapper<3, random::gpu::rand48_rng>;
template class brownian_wrapper<2, random::gpu::mrg32k3a_rng>;
template class brownian_wrapper<3, random::gpu::mrg32k3a_rng>;

} // namespace integrators
} // namespace gpu
} // namespace mdsim
} // namespace halmd
