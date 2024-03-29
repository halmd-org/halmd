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

#ifndef HALMD_RANDOM_GPU_RANDOM_KERNEL_HPP
#define HALMD_RANDOM_GPU_RANDOM_KERNEL_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>

namespace halmd {
namespace random {
namespace gpu {

template <typename RandomNumberGenerator>
struct random_wrapper
{
    cuda::function<void (float*, unsigned int, RandomNumberGenerator)> uniform;
    cuda::function<void (unsigned int*, unsigned int, RandomNumberGenerator)> get;
    cuda::function<void (float*, unsigned int, float, float, RandomNumberGenerator)> normal;
    static random_wrapper kernel;
};

template <typename RandomNumberGenerator>
random_wrapper<RandomNumberGenerator>& get_random_kernel()
{
    return random_wrapper<RandomNumberGenerator>::kernel;
}

} // namespace random
} // namespace gpu
} // namespace halmd

#endif /* ! HALMD_RANDOM_GPU_RANDOM_KERNEL_HPP */
