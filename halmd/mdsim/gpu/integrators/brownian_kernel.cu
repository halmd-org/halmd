/*
 * Copyright © 2016      Manuel Dibak
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

#include <cmath>

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
static texture<float4> param_;

template <int dimension, typename float_type, typename vector_type>
__device__ void update_displacement(
    float_type const diff_const_par
  , float_type const diff_const_perp
  , float_type const prop_strength
  , vector_type& r
  , vector_type const& u
  , vector_type const& f
  , float timestep
  , float temperature
  , float_type eta1
  , float_type eta2
  , float_type eta3
)
{
    vector_type dr_r, dr_s;

    // random part
    dr_r[0] = eta1;
    dr_r[1] = eta2;
    if(dimension % 2) {
        dr_r[2] = eta3;
    }

    // systeatic part
    float_type f_u = inner_prod(f, u);
    dr_s = (timestep * f_u * diff_const_par / temperature + prop_strength * timestep) * u +
        (timestep * diff_const_perp / temperature) * (f - f_u * u);

    r += dr_r + dr_s;
}

template <typename float_type>
__device__ void update_orientation(
    float_type const diff_const_rot
  , fixed_vector<float_type, 2>& u
  , fixed_vector<float_type, 2> const& tau
  , float timestep
  , float temperature
  , float_type eta
  , float_type ignore
)
{
    fixed_vector<float_type, 2> e;

    // construct trihedron along particle orientation for movement in x-y plane
    // e1 lies in x-y plane
    e[0] = -u[1];
    e[1] = u[0];

    // first term is the random torque, second is the systematic torque
    float_type omega = eta + tau[0] * diff_const_rot * timestep / temperature;

    u = __cosf(omega) * u + __sinf(omega) * e;

    // ensure normalization, TODO: is this really necessary?
    if ((float) norm_2(u) > 2e-38){
        u /= norm_2(u);
    }
}

template <typename float_type>
__device__ void update_orientation(
    float_type const diff_const_rot
  , fixed_vector<float_type, 3>& u
  , fixed_vector<float_type, 3> const& tau
  , float timestep
  , float temp
  , float_type eta1
  , float_type eta2
)
{
    fixed_vector<float_type, 3> e1, e2;

    //construct trihedron along particle orientation
    if ((float) u[1] > 2e-38 || (float) u[2] > 2e-38) {
        e1[0] = 0;
        e1[1] = u[2];
        e1[2] = -u[1];
    } else {
        e1[0] = u[1];
        e1[1] = -u[0];
        e1[2] = 0;
    }
    e1 /= norm_2(e1);
    e2 = cross_prod(u, e1);
    e2 /= norm_2(e2);

    // first two terms are the random angular velocity, the lst part is the systematic torque
    fixed_vector<float_type, 3> omega = eta1 * e1 + eta2 * e2 + tau * diff_const_rot * timestep / temp ;
    float omega_abs = norm_2(omega);
    omega /= omega_abs;
    // Ω = eta1 * e1 + eta2 * e2
    // => Ω × u = (eta1 * e1 × u + eta2  * e2 × u) = eta2 * e1 - eta1 * e2
    u = (1 - __cosf(omega_abs)) * inner_prod(omega, u) * omega +
        __cosf(omega_abs) * u +
        __sinf(omega_abs) * cross_prod(omega, u);

    // ensure normalization
    if ((float) norm_2(u) != 0){
        u /= norm_2(u);
    }
}

template <
    int dimension
  , typename float_type
  , typename ptr_type
  , typename const_ptr_type
  , typename gpu_vector_type
  , typename gpu_pseudo_vector_type
  , typename rng_type
>
__global__ void integrate(
    ptr_type g_position
  , ptr_type g_orientation
  , gpu_vector_type* g_image
  , const_ptr_type g_velocity
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

    // read position, species, velocity, mass, image from global memory
    vector_type r, u, v;
    unsigned int species;
    float mass;

    float_type rng_disp_cache = 0;
    bool rng_disp_cached = false;
    float_type rng_rot_cache = 0;
    bool rng_rot_cached = false;

    //read random number generator state from global device memory
    typename rng_type::state_type state = rng[thread];

    for (unsigned int i = thread; i < nparticle; i += nthread) {
        // read in relevant variables
        tie(r, species) <<= g_position[i];
        tie(u, mass) <<= g_orientation[i];

        // normal distribution parameters
        fixed_vector<float, 4> diff_const = tex1Dfetch(param_, species);
        float_type const diff_const_perp  = diff_const[0];
        float_type const diff_const_par   = diff_const[1];
        float_type const diff_const_rot   = diff_const[2];
        float_type const prop_str   = diff_const[3];
        float_type const sigma_perp = sqrtf(2 * diff_const_perp * timestep );
        float_type const sigma_par  = sqrtf(2 * diff_const_par * timestep );
        float_type const sigma_rot  = sqrtf(2 * diff_const_rot * timestep );

        vector_type f   = static_cast<float_vector_type>(g_force[i]);
        // TODO: torque type
        vector_type tau = static_cast<float_vector_type>(g_torque[i]);

        // draw random numbers
        float_type eta1, eta2, eta3;
        tie(eta1, eta2) =  random::gpu::normal(rng, state, 0, sigma_perp);
        if (dimension % 2 == 1) {
            if (rng_disp_cached) {
                eta3 = rng_disp_cache;
            } else {
                tie(eta3, rng_disp_cache) = random::gpu::normal(rng, state, 0, sigma_par);
            }
            rng_disp_cached = !rng_disp_cached;
        }

        // Brownian integration
        update_displacement<dimension, float_type, vector_type>(
            diff_const_par, diff_const_perp, prop_str, r, u, f, timestep, temp, eta1, eta2, eta3
        );

        // enforce periodic boundary conditions
        float_vector_type image = box_kernel::reduce_periodic(r, box_length);

        if (dimension == 2) {
            if (rng_rot_cached) {
                eta1 = rng_rot_cache;
            } else {
                tie(eta1, rng_rot_cache) =  random::gpu::normal(rng, state, 0, sigma_rot);
            }
            rng_rot_cached = !rng_rot_cached;
        } else {
            tie(eta1, eta2) =  random::gpu::normal(rng, state, 0, sigma_rot);
        }

        // update orientation last (Ito interpretation)
        update_orientation(diff_const_rot, u, tau, timestep, temp, eta1, eta2);

        // store position, species, image in global memory
        g_position[i] <<= tie(r, species);
        g_orientation[i] <<= tie(u, mass);

        if (!(image == float_vector_type(0))) {
            g_image[i] = image + static_cast<float_vector_type>(g_image[i]);
        }
    }
    // store random number generator state in global device memory
    rng[thread] = state;
}

} // namespace brownian_kernel

template <int dimension, typename float_type, typename rng_type>
cuda::texture<float4> brownian_wrapper<dimension, float_type, rng_type>::param = brownian_kernel::param_;

template <int dimension, typename float_type, typename rng_type>
brownian_wrapper<dimension, float_type, rng_type> const
brownian_wrapper<dimension, float_type, rng_type>::kernel = {
    brownian_kernel::integrate<dimension, float_type, ptr_type, const_ptr_type>
};

// explicit instantiation
#ifdef USE_GPU_SINGLE_PRECISION
template class brownian_wrapper<2, float, random::gpu::rand48_rng>;
template class brownian_wrapper<3, float, random::gpu::rand48_rng>;
template class brownian_wrapper<2, float, random::gpu::mrg32k3a_rng>;
template class brownian_wrapper<3, float, random::gpu::mrg32k3a_rng>;
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
template class brownian_wrapper<2, dsfloat, random::gpu::rand48_rng>;
template class brownian_wrapper<3, dsfloat, random::gpu::rand48_rng>;
template class brownian_wrapper<2, dsfloat, random::gpu::mrg32k3a_rng>;
template class brownian_wrapper<3, dsfloat, random::gpu::mrg32k3a_rng>;
#endif

} // namespace integrators
} // namespace gpu
} // namespace mdsim
} // namespace halmd
