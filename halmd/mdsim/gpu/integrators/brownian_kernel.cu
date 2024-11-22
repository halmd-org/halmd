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

/** array of diffusion constants */
static texture<float2> param_;

template <
    int dimension
  , typename float_type
  , typename ptr_type
  , typename gpu_vector_type
  , typename rng_type
>
__global__ void integrate(
    ptr_type g_position
  , gpu_vector_type* g_image
  , gpu_vector_type const* g_force
  , float timestep
  , float temp
  , rng_type rng
  , unsigned int nparticle
  , fixed_vector<float, dimension> box_length
)
{
    typedef fixed_vector<float_type, dimension> vector_type;
    typedef fixed_vector<float, dimension> float_vector_type;

    // kernel execution parameters
    unsigned int const thread = GTID;
    unsigned int const nthread = GTDIM;

    // read position, species, image from global memory
    vector_type r;
    unsigned int species;
    float mass;

    float_type rng_disp_cache = 0;
    bool rng_disp_cached = false;

    //read random number generator state from global device memory
    typename rng_type::state_type state = rng[thread];

    for (unsigned int i = thread; i < nparticle; i += nthread) {
        // read position (do this either way because we need the species)
        tie(r, species) <<= g_position[i];

        // diffusion constants (we will need at least one of those)
        fixed_vector<float, 2> diff_const = tex1Dfetch(param_, species);

        vector_type f = static_cast<float_vector_type>(g_force[i]);

        float_type const diff_const_disp = diff_const[0];
        float_type const sigma_disp = sqrtf(2 * diff_const_disp * timestep);

        // draw Gaussian random vector
        vector_type dr;
        tie(dr[0], dr[1]) =  random::gpu::normal(rng, state, 0, sigma_disp);
        if (dimension == 3) {
            if (rng_disp_cached) {
                dr[2] = rng_disp_cache;
            } else {
                tie(dr[2], rng_disp_cache) = random::gpu::normal(rng, state, 0, sigma_disp);
            }
            rng_disp_cached = !rng_disp_cached;
        }

        // Brownian integration: Euler-Maruyama scheme
        r += dr + (diff_const * timestep / temperature) * f;

        // enforce periodic boundary conditions
        float_vector_type image = box_kernel::reduce_periodic(r, box_length);

        // store position and image (do this here because the orientation doesn't change the position or species)
        g_position[i] <<= tie(r, species);
        if (!(image == float_vector_type(0))) {
            g_image[i] = image + static_cast<float_vector_type>(g_image[i]);
        }
    }

    // store random number generator state in global device memory
    rng[thread] = state;
}

} // namespace brownian_kernel

template <int dimension, typename float_type, typename rng_type>
cuda::texture<float2> brownian_wrapper<dimension, float_type, rng_type>::param = brownian_kernel::param_;

template <int dimension, typename float_type, typename rng_type>
brownian_wrapper<dimension, float_type, rng_type> const
brownian_wrapper<dimension, float_type, rng_type>::kernel = {
    brownian_kernel::integrate<dimension, float_type, ptr_type>
};

// explicit instantiation
#ifdef USE_GPU_SINGLE_PRECISION
template class brownian_wrapper<2, float, random::gpu::rand48_rng>;
template class brownian_wrapper<3, float, random::gpu::rand48_rng>;
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
template class brownian_wrapper<2, dsfloat, random::gpu::rand48_rng>;
template class brownian_wrapper<3, dsfloat, random::gpu::rand48_rng>;
#endif

} // namespace integrators
} // namespace gpu
} // namespace mdsim
} // namespace halmd
