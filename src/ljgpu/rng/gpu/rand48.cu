/* Parallelized rand48 random number generator for CUDA
 *
 * Copyright © 2007-2009  Peter Colberg
 *
 * This program is free software: you can redistribute it and/or modify
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

#include <ljgpu/algorithm/gpu/base.cuh>
#include <ljgpu/rng/gpu/rand48.cuh>
#include <ljgpu/rng/gpu/rand48.hpp>
using namespace ljgpu::gpu::rand48;

namespace ljgpu { namespace cu { namespace rand48
{

/*
 * This is a parallel version of the Unix rand48 generator for CUDA.
 * It is based on the GNU Scientific Library rand48 implementation.
 */

/**
 * compute leapfrog multipliers for initialization
 */
__global__ void leapfrog(uint48* g_la)
{
    const uint48 a(0xE66D, 0xDEEC, 0x0005);

    //
    // leapfrog multiplier:
    //   A = a^N mod m
    //

    uint48 A = a;
    uint48 x = a;

    // fast exponentiation by squares
    for (uint k = GTID; k > 0; k >>= 1) {
        if (k % 2 == 1) {
            A = muladd(x, A, 0);
        }
        x = muladd(x, x, 0);
    }

    g_la[GTID] = A;
}

/**
 * initialize generator with 32-bit integer seed
 */
__global__ void set(uint48 const* g_la, uint48 const* g_lc, uint48 *g_a, uint48 *g_c, uint seed)
{
    const uint48 c(0x000B, 0, 0);

    //
    // leapfrog addend:
    //   C = (c * sum(n = 0..(N-1), a^n)) mod m
    //

    const uint48 A = g_la[GTID];
    const uint48 C = muladd(c, g_lc[GTID], c);

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

/**
 * restore generate state
 */
__global__ void restore(uint48 const* g_la, uint48 const* g_lc, uint48 *g_a, uint48 *g_c, ushort3 state)
{
    const uint48 c(0x000B, 0, 0);

    const uint48 A = g_la[GTID];
    const uint48 C = muladd(c, g_lc[GTID], c);

    if (GTID == GTDIM - 1) {
        // store leapfrog constants
        *g_a = A;
        *g_c = C;

        g_state[0] = state;
    }
    else {
        // generate initial states
        g_state[GTID + 1] = muladd(A, state, C);
    }
}

/**
 * save generator state
 */
__global__ void save(ushort3 *state)
{
    if (GTID == 0) {
        *state = g_state[0];
    }
}

/**
 * fill array with uniform random numbers in [0.0, 1.0)
 */
__global__ void uniform(float* v, uint len)
{
    ushort3 x = g_state[GTID];

    for (uint k = GTID; k < len; k += GTDIM) {
        v[k] = uniform(x);
    }

    g_state[GTID] = x;
}

/**
 * returns random integer in [0, 2^32-1]
 */
__device__ uint get(ushort3& state)
{
    uint r = (state.z << 16UL) + state.y;
    state = muladd(a, state, c);
    return r;
}

/**
 * fill array with random integers in [0, 2^32-1]
 */
__global__ void get(uint* v, uint len)
{
    ushort3 x = g_state[GTID];

    for (uint k = GTID; k < len; k += GTDIM) {
        v[k] = get(x);
    }

    g_state[GTID] = x;
}

}}} // namespace ljgpu::cu::rand48

namespace ljgpu { namespace gpu
{

/**
 * device function wrappers
 */
cuda::function<void (uint48*)>
    rand48::leapfrog(cu::rand48::leapfrog);
cuda::function<void (uint48 const*, uint48 const*, uint48*, uint48*, uint)>
    rand48::set(cu::rand48::set);
cuda::function<void (uint48 const*, uint48 const*, uint48*, uint48*, ushort3)>
    rand48::restore(cu::rand48::restore);
cuda::function<void (ushort3*)>
    rand48::save(cu::rand48::save);
cuda::function<void (float*, uint)>
    rand48::uniform(cu::rand48::uniform);
cuda::function<void (uint*, uint)>
    rand48::get(cu::rand48::get);

/**
 * device constant wrappers
 */
cuda::symbol<uint48> rand48::a(cu::rand48::a);
cuda::symbol<uint48> rand48::c(cu::rand48::c);
cuda::symbol<ushort3*> rand48::state(cu::rand48::g_state);

}} // namespace ljgpu::gpu
