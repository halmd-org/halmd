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

#ifndef HALMD_RANDOM_GPU_MRG32K3A_KERNEL_HPP
#define HALMD_RANDOM_GPU_MRG32K3A_KERNEL_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>

#include <halmd/random/gpu/mrg32k3a_kernel.cuh>

namespace halmd {
namespace random {
namespace gpu {

/**
 * CUDA C++ wrapper
 */
struct mrg32k3a_wrapper
{
    typedef mrg32k3a_rng rng_type;

    cuda::function<void (rng_type::state_type*, unsigned int)> seed;
    static mrg32k3a_wrapper kernel;
};

// syntactic sugar
inline mrg32k3a_wrapper& get_mrg32k3a_kernel()
{
    return mrg32k3a_wrapper::kernel;
}

} // namespace random
} // namespace gpu
} // namespace halmd

#endif /* ! HALMD_RANDOM_GPU_MRG32K3A_KERNEL_HPP */
