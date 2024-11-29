/*
 * Copyright © 2024 Felix Höfling
 * Copyright © 2024 Max Orteu
 * Copyright © 2016 Manuel Dibak
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

template <
    int dimension
  , typename float_type
  , typename ptr_type
  , typename gpu_vector_type
  , typename rng_type
>
__global__ void integrate(
    cudaTextureObject_t t_param
  , ptr_type g_position
  , gpu_vector_type* g_image
  , gpu_vector_type const* g_force
  , float timestep
  , float temperature
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

    float_type rng_cache = 0;
    bool rng_cached = false;

    //read random number generator state from global device memory
    typename rng_type::state_type state = rng[thread];

    for (unsigned int i = thread; i < nparticle; i += nthread) {
        // read position and species from global memory
        vector_type r;
        unsigned int species;
        tie(r, species) <<= g_position[i];

        // read force from global memory
        vector_type f = static_cast<float_vector_type>(g_force[i]);

        // read diffusion constant from texture
        float diffusion = tex1Dfetch<float>(t_param, species);
        float_type const sigma = sqrtf(2 * diffusion * timestep);

        // draw Gaussian random vector
        vector_type dr;
        tie(dr[0], dr[1]) =  random::gpu::normal(rng, state, 0, sigma);
        if (dimension == 3) {
            if (rng_cached) {
                dr[2] = rng_cache;
            } else {
                tie(dr[2], rng_cache) = random::gpu::normal(rng, state, 0, sigma);
            }
            rng_cached = !rng_cached;
        }

        // Brownian integration: Euler-Maruyama scheme
        r += dr + (diffusion * timestep / temperature) * f;

        // enforce periodic boundary conditions
        float_vector_type image = box_kernel::reduce_periodic(r, box_length);

        // store position and image in global memory
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
brownian_wrapper<dimension, float_type, rng_type>
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
