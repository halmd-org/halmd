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
#include <cmath>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace integrators {
namespace brownian_kernel {

/**array of diffusion constants */
static texture<float4> param_;

/**cross product of two vectors

template <typename vector_type>
__device__ vector_type cross_prod(vector_type a, vector_type b)
{
    vector_type c;
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
    return c;
}

/** inner product of two vectors

template <typename float_type, typename vector_type>
__device__ float_type inner_prod(vector_type a, vector_type b)
{
    return a[0] * b[0] 
}

#ifdef USE_VERLET_DSFUN
template __device__ fixed_vector<dsfloat, 3> cross_prod<3, fixed_vector<dsfloat, 3> >(fixed_vector<dsfloat, 3>* a, fixed_vector<dsfloat, 3>* b);
#else
template __device__ fixed_vector<float, 3> cross_prod<3, fixed_vector<float, 3> >(fixed_vector<float, 3>* a, fixed_vector<float, 3>* b);
#endif

template <typename float_type, typename vector_type>
__device__ float_type norm_2(vector_type* a)
{
    return a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
}
*/
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
  , float4* g_orientation
  , gpu_vector_type* g_image
  , float4 const* g_velocity
  , gpu_vector_type const* g_force
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
    float const mean = 0.f;
    
    //read random number generator state from global device memory
    typename rng_type::state_type state = rng[thread];

    for (uint i = thread; i < nparticle; i += nthread) {

#ifdef USE_VERLET_DSFUN
        tie(r, species) <<= tie(g_position[i], g_position[i + nplace]);
        tie(u, nothing) <<= tie(g_orientation[i], g_orientation[i + nplace]);
        tie(v, mass) <<= tie(g_velocity[i], g_velocity[i + nplace]);
#else
        tie(r, species) <<= g_position[i];
        tie(u, nothing) <<= g_orientation[i];
        tie(v, mass) <<= g_velocity[i];
#endif
        //initialize random displacement and trihedron
        vector_type dr, e1, e2;
        
        //construct trihedron along particle orientation
        if ( (float) u[1] > epsilon || (float) u[2] > epsilon) {
            e1[0] = 0; e1[1] = u[2]; e1[2] = -u[1];
        }
        else {
            e1[0] = u[1]; e1[1] = -u[0]; e1[2] = 0;
        }
        e1 /= norm_2(e1);
        // compute cross product u x e1
        e2 = cross_prod(u, e1);
        e2 /= norm_2(e2);
    
        //normal distribution parameters
        fixed_vector<float, 4> D = tex1Dfetch(param_, species);
        float_type const D_perp =  D[0];
        float_type const D_par = D[1];
        float_type const D_rot = D[2];
        float_type const prop_str = D[3];
        float_type const sigma_perp = sqrtf(2 * D_perp * timestep);
        float_type const sigma_par = sqrtf(2 * D_par * timestep);
        float_type const sigma_rot = sqrtf(2 * D_rot * timestep);
   
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
        // random displacement
        dr = eta1 * e1 + eta2 * e2 + (eta3 + prop_str * timestep) * u; 
        vector_type f = static_cast<float_vector_type>(g_force[i]);
        //systematic displacement
        float_type f_u = inner_prod(f, u);
        //r += (f_u * D_par / temp + prop_str) * u + (D_perp / temp) * (f - f_u * u);
        r += dr; 
        // update orientation last (Ito interpretation)
        // no torque included!
        tie(eta1, eta2) =  random::gpu::normal(rng, state, mean, sigma_rot);
        vector_type omega = eta1 * e1 + eta2 * e2;
        float const abs = norm_2(omega);
        // Ω = eta1 * e1 + eta2 * e2
        // => Ω × u = (eta1 * e1 × u + eta2  * e2 × u) = eta2 * e1 - eta1 * e2
        u = cos(abs) * u + sin(abs) / abs * cross_prod(u, omega);
        //ensure normalization
        u /= norm_2(u);
        // enforce periodic boundary conditions
        float_vector_type image = box_kernel::reduce_periodic(r, box_length);

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
template class brownian_wrapper<3, random::gpu::rand48_rng>;
template class brownian_wrapper<3, random::gpu::mrg32k3a_rng>;
// let's stay 3 dimensional for now
//template class brownian_wrapper<2, random::gpu::rand48_rng>;

} // namespace integrators
} // namespace gpu
} // namespace mdsim
} // namespace halmd
