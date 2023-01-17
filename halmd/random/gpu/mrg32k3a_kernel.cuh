/*
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

#ifndef HALMD_RANDOM_GPU_MRG32K3A_KERNEL_CUH
#define HALMD_RANDOM_GPU_MRG32K3A_KERNEL_CUH

#include <curand_kernel.h>

#include <halmd/config.hpp>

namespace halmd {
namespace random {
namespace gpu {

struct mrg32k3a_rng
{
    /** per-thread generator state type */
    typedef curandStateMRG32k3a state_type;

    HALMD_GPU_ENABLED state_type& operator[](unsigned int thread) const
    {
        return g_state[thread];
    }
    /** generator states in global device memory */
    state_type* g_state;
};

/**
 * returns random integer in [0, 2^32-1]
 */
inline __device__ unsigned int get(mrg32k3a_rng const& rng, mrg32k3a_rng::state_type& state)
{
    unsigned int variate = curand(&state);
    return variate;
}

} // namespace random
} // namespace gpu
} // namespace halmd

#endif /* ! HALMD_RANDOM_GPU_MRG32K3A_KERNEL_CUH */
