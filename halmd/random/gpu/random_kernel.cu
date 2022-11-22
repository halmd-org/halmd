/*
 * Copyright © 2010  Peter Colberg
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

#include <halmd/algorithm/gpu/scan_kernel.cuh>
#include <halmd/random/gpu/normal_distribution.cuh>
#include <halmd/random/gpu/random_number_generator.cuh>
#include <halmd/random/gpu/random_kernel.hpp>
#include <halmd/utility/gpu/thread.cuh>

namespace halmd {
namespace random {
namespace gpu {
namespace random_kernel {

/**
 * fill array with uniform random numbers in [0.0, 1.0)
 */
template <typename Rng>
__global__ void uniform(float* v, unsigned int len, Rng rng)
{
    typename Rng::state_type state = rng[GTID];

    for (unsigned int k = GTID; k < len; k += GTDIM) {
        v[k] = uniform(rng, state);
    }

    rng[GTID] = state;
}

/**
 * fill array with random integers in [0, 2^32-1]
 */
template <typename Rng>
__global__ void get(unsigned int* v, unsigned int len, Rng rng)
{
    typename Rng::state_type state = rng[GTID];

    for (unsigned int k = GTID; k < len; k += GTDIM) {
        v[k] = get(rng, state);
    }

    rng[GTID] = state;
}

/**
 * fill array with normal distributed random numbers in [0.0, 1.0)
 */
template <typename Rng>
__global__ void normal(float* v, unsigned int len, float mean, float sigma, Rng rng)
{
    typename Rng::state_type state = rng[GTID];

    for (unsigned int k = 2 * GTID; k < len - 1; k += 2 * GTDIM) {
        tie(v[k], v[k + 1]) = normal(rng, state, mean, sigma);
    }

    // if len is odd we need to still write the last element (in the first thread)
    if (GTID == 0 && len % 2 == 1) {
        float discard;
        tie(v[len - 1], discard) = normal(rng, state, mean, sigma);
    }

    rng[GTID] = state;
}


} // namespace random_kernel

/**
 * CUDA C++ wrappers
 */
template <typename Rng>
random_wrapper<Rng> random_wrapper<Rng>::kernel = {
    random_kernel::uniform<Rng>
  , random_kernel::get<Rng>
  , random_kernel::normal<Rng>
};

template class random_wrapper<random::gpu::rand48_rng>;

} // namespace random
} // namespace gpu
} // namespace halmd
