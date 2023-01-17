/*
 * Copyright © 2010  Peter Colberg
 * Copyright © 2023  Jaslo Ziska
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

#ifndef HALMD_RANDOM_GPU_RANDOM_KERNEL_CUH
#define HALMD_RANDOM_GPU_RANDOM_KERNEL_CUH

#include <halmd/random/gpu/rand48_kernel.cuh>
#include <halmd/random/gpu/mrg32k3a_kernel.cuh>

/**
 * Returns uniform random number in [0.0, 1.0)
 *
 * This function generates uniform and unbiased floats in [0.0, 1.0) from the 24 most significant bits of
 * an unsigned integers returned by the RNG using the procedure described in this blog post by Daniel Lemire:
 * <https://lemire.me/blog/2017/02/28/how-many-floating-point-numbers-are-in-the-interval-01/>
 */
template <typename rng_type>
inline __device__ float uniform(rng_type const& rng, typename rng_type::state_type& state)
{
    unsigned int variate = get(rng, state);
    float result = static_cast<float>(variate >> 8) / 16777216.f;
    return result;
}

#endif /* ! HALMD_RANDOM_GPU_RANDOM_KERNEL_CUH */
