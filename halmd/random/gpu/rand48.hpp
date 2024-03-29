/*
 * Copyright © 2007-2010  Peter Colberg
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

#ifndef HALMD_RANDOM_GPU_RAND48_HPP
#define HALMD_RANDOM_GPU_RAND48_HPP

#include <algorithm>
#include <iostream>

#include <halmd/algorithm/gpu/scan.hpp>
#include <halmd/random/gpu/rand48_kernel.hpp>

namespace halmd {
namespace random {
namespace gpu {

/**
 * Parallelized rand48 random number generator for CUDA
 */
class rand48
{
public:
    typedef rand48_rng rng_type;

    static char const* name() {
        return "rand48";
    }

    /**
     * initialize random number generator with CUDA execution dimensions
     */
    rand48(dim3 blocks, dim3 threads)
      : dim(blocks, threads)
      , g_state_(dim.threads())
    {}

    /**
     * seed generator with 32-bit integer
     */
    void seed(unsigned int value)
    {
        // compute leapfrog multipliers for initialization
        cuda::memory::device::vector<uint48> g_A(dim.threads()), g_C(dim.threads());
        get_rand48_kernel().leapfrog.configure(dim.grid, dim.block);
        get_rand48_kernel().leapfrog(g_A);

        // compute leapfrog addends for initialization
        cuda::copy(g_A.begin(), g_A.end(), g_C.begin());
        algorithm::gpu::scan<uint48> scan(g_C.size(), dim.threads_per_block());
        scan(g_C);

        // initialize generator with seed
        cuda::memory::device::vector<uint48> g_a(1), g_c(1);
        cuda::memory::host::vector<uint48> h_a(1), h_c(1);
        get_rand48_kernel().seed.configure(dim.grid, dim.block);
        get_rand48_kernel().seed(g_A, g_C, g_a, g_c, g_state_, value);
        cuda::copy(g_a.begin(), g_a.end(), h_a.begin());
        cuda::copy(g_c.begin(), g_c.end(), h_c.begin());

        // set leapfrog constants for constant device memory
        rng_.a = h_a.front();
        rng_.c = h_c.front();
        rng_.g_state = g_state_.data();
    }

    cuda::config const dim;

    rand48_rng const& rng() const
    {
        return rng_;
    }

private:
    cuda::memory::device::vector<ushort3> g_state_;
    rand48_rng rng_;
};

} // namespace random
} // namespace gpu
} // namespace halmd

#endif /* ! HALMD_RANDOM_GPU_RAND48_HPP */
