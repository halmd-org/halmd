/*
 * Copyright © 2010  Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef HALMD_RANDOM_GPU_RANDOM_KERNEL_CUH
#define HALMD_RANDOM_GPU_RANDOM_KERNEL_CUH

#include <halmd/random/gpu/rand48_kernel.cuh>

namespace halmd
{
namespace random { namespace gpu
{

union random_number_generator
{
private:
    /**
     * Available random number generators.
     */
    rand48_rng rand48_;
    // gfsr4_rng gfsr4_;

public:
    /**
     * These constructors are used on the host to cuda::copy the
     * parameters of a specific random number generator to the GPU.
     */
    __host__ random_number_generator(rand48_rng const& rng) : rand48_(rng) {}
    // __host__ random_number_generator(gfsr4_rng const& rng) : gfsr4_(rng) {}

    /**
     * Default constructor for __constant__ definition.
     */
    __device__ random_number_generator() {}
};

/**
 * Select random number generator in __device__ function.
 */
template <typename RandomNumberGenerator>
__device__ __host__ RandomNumberGenerator const& get(random_number_generator const& rng)
{
    return reinterpret_cast<RandomNumberGenerator const&>(rng);
}

}} // namespace random::gpu

} // namespace halmd

#endif /* ! HALMD_RANDOM_GPU_RANDOM_KERNEL_CUH */
