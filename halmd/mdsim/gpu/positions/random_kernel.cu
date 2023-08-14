/*
 * Copyright © 2017 Arthur Straube
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

#include <halmd/mdsim/gpu/positions/random_kernel.hpp>
#include <halmd/random/gpu/random_number_generator.cuh>
#include <halmd/utility/gpu/thread.cuh>

// using namespace halmd::algorithm::gpu;

//
// Random (uniform) distribution of particles over the simulation box
//

namespace halmd {
namespace mdsim {
namespace gpu {
namespace positions {
namespace random_kernel {

/**
 * generate random positions within slab centred around the origin
 *
 * work in single precision, reset high precision part of dsfloat
 */
template <typename vector_type, typename rng_type>
__global__ void uniform(
    float4* g_r
  , unsigned int npart
  , unsigned int nplace
  , vector_type slab_length
  , rng_type rng
)
{
    enum { dimension = vector_type::static_size };

    // read random number generator state from global device memory
    typename rng_type::state_type state = rng[GTID];

    for (uint i = GTID; i < npart; i += GTDIM) {
        vector_type r;
        unsigned int species;

        tie(r, species) <<= g_r[i];

        for (uint j = 0; j < dimension; ++j) {
            r[j] = uniform(rng, state) - .5f;
        }
        r = element_prod(r, slab_length);

        g_r[i] <<= tie(r, species);
#ifdef USE_VERLET_DSFUN
        g_r[i + nplace] = fixed_vector<float, 4>(0);
#endif
    }

    // store random number generator state in global device memory
    rng[GTID] = state;
}

} // namespace random_kernel

template <int dimension, typename rng_type>
random_wrapper<dimension, rng_type> random_wrapper<dimension, rng_type>::kernel = {
    random_kernel::uniform<fixed_vector<float, dimension>, rng_type>
};

template class random_wrapper<3, random::gpu::rand48_rng>;
template class random_wrapper<2, random::gpu::rand48_rng>;

} // namespace positions
} // namespace gpu
} // namespace mdsim
} // namespace halmd
