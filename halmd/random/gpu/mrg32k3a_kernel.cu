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

#include <halmd/random/gpu/mrg32k3a_kernel.hpp>
#include <halmd/utility/gpu/thread.cuh>

namespace halmd {
namespace random {
namespace gpu {
namespace mrg32k3a_kernel {

__global__ void seed(curandStateMRG32k3a* state, uint seed)
{
    unsigned int id = GTID;
    /* Each thread gets same seed, a different sequence
       number, no offset */
    curand_init(seed, id, 0, &state[id]);
}

} // namespace mrg32k3a_kernel

/**
 * CUDA C++ wrappers
 */
mrg32k3a_wrapper mrg32k3a_wrapper::kernel = {
  mrg32k3a_kernel::seed
};

} // namespace random
} // namespace gpu
} // namespace halmd
