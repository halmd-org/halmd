/*
 * Copyright © 2007-2010  Peter Colberg
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

//
// This is a parallel version of the Unix rand48 generator for CUDA.
// It is based on the rand48 generator of the GNU Scientific Library.
// The file rng/rand48.c was written by James Theiler and Brian Gough
// and is licensed under the GPL v3 or later.
//

#include <halmd/algorithm/gpu/scan_kernel.cuh>
#include <halmd/random/gpu/rand48_kernel.hpp>
#include <halmd/utility/gpu/thread.cuh>

namespace halmd {
namespace random {
namespace gpu {
namespace rand48_kernel {

/*
 * This is a parallel version of the Unix rand48 generator for CUDA.
 * It is based on the GNU Scientific Library rand48 implementation.
 */

/**
 * compute leapfrog multipliers for initialization
 */
__global__ void leapfrog(uint48* g_A)
{
    uint48 const a = make_uint48(0xE66D, 0xDEEC, 0x0005);
    uint48 const zero = make_uint48(0, 0, 0);

    //
    // leapfrog multiplier:
    //   A = a^N mod m
    //

    uint48 A = a;
    uint48 x = a;

    // fast exponentiation by squares
    for (uint k = GTID; k > 0; k >>= 1) {
        if (k % 2 == 1) {
            A = muladd(x, A, zero);
        }
        x = muladd(x, x, zero);
    }

    g_A[GTID] = A;
}

/**
 * initialize generator with 32-bit integer seed
 */
__global__ void seed(uint48 const* g_A, uint48 const* g_C, uint48 *g_a, uint48 *g_c, ushort3* g_state, uint seed)
{
    uint48 const c = make_uint48(0x000B, 0, 0);

    //
    // leapfrog addend:
    //   C = (c * sum(n = 0..(N-1), a^n)) mod m
    //

    uint48 const A = g_A[GTID];
    uint48 const C = muladd(c, g_C[GTID], c);

    if (GTID == GTDIM - 1) {
        // store leapfrog constants
        *g_a = A;
        *g_c = C;
    }

    // default seed
    ushort3 x = make_ushort3(0x330E, 0xABCD, 0x1234);

    if (seed > 0) {
        x.y = seed & 0xFFFF;
        x.z = (seed >> 16) & 0xFFFF;
    }

    // generate initial state
    g_state[GTID] = muladd(A, x, C);
}

} // namespace rand48_kernel

/**
 * CUDA C++ wrappers
 */
rand48_wrapper const rand48_wrapper::kernel = {
    rand48_kernel::leapfrog
  , rand48_kernel::seed
};

} // namespace random
} // namespace gpu

template class algorithm::gpu::scan_wrapper<uint48>;

} // namespace halmd
