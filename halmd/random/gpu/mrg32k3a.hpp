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

#ifndef HALMD_RANDOM_GPU_MRG32K3A_HPP
#define HALMD_RANDOM_GPU_MRG32K3A_HPP

#include <algorithm>
#include <iostream>

#include <halmd/random/gpu/mrg32k3a_kernel.hpp>

namespace halmd {
namespace random {
namespace gpu {

/**
 * Implement the MRG32K3a random number generator from the cuRAND library
 */
class mrg32k3a
{
public:
    typedef mrg32k3a_rng rng_type;

    static char const* name() {
        return "MRG32k3a";
    }

    /**
     * initialize random number generator with CUDA execution dimensions
     */
    mrg32k3a(dim3 blocks, dim3 threads)
      : dim(blocks, threads)
      , g_state_(dim.threads())
    {}

    /**
     * setup kernel with seed 
     */
    void seed(unsigned int value)
    {
        cuda::configure(dim.grid, dim.block);
        // initialize generator with seed
        get_mrg32k3a_kernel().seed(g_state_, value);

        rng_.g_state = g_state_.data();
    }

    cuda::config const dim;

    mrg32k3a_rng const& rng() const
    {
        return rng_;
    }

private:
    cuda::vector<curandStateMRG32k3a> g_state_;
    mrg32k3a_rng rng_;
};

} // namespace random
} // namespace gpu
} // namespace halmd

#endif /* ! HALMD_RANDOM_GPU_MRG32K3A_HPP */
