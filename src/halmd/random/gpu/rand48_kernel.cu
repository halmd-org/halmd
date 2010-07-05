/*
 * Copyright © 2007-2010  Peter Colberg and Felix Höfling
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

#include <halmd/random/gpu/rand48_kernel.cuh>
#include <halmd/utility/gpu/thread.cuh>

namespace halmd
{
namespace random { namespace gpu
{

//
// This is a parallel version of the Unix rand48 generator for CUDA.
// It is based on the GNU Scientific Library rand48 implementation.
//

typedef ushort3 state_type;

/** leapfrogging multiplier */
__constant__ uint48 a;
/** leapfrogging addend */
__constant__ uint48 c;
/** generator state in global device memory */
__constant__ state_type* g_state;


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
 * returns uniform random number in [0.0, 1.0)
 */
__device__ float uniform(state_type& state)
{
    float r = state.z / 65536.f + state.y / 4294967296.f;
    state = muladd(a, state, c);
    return r;
}

/**
 * generate 2 random numbers from Gaussian distribution with given variance
 */
__device__ void normal(float& r1, float& r2, float var, state_type& state)
{
    //
    // The Box-Muller transformation for generating random numbers
    // in the normal distribution was originally described in
    //
    // G.E.P. Box and M.E. Muller, A Note on the Generation of
    // Random Normal Deviates, The Annals of Mathematical Statistics,
    // 1958, 29, p. 610-611
    //
    // Here, we use instead the faster polar method of the Box-Muller
    // transformation, see
    //
    // D.E. Knuth, Art of Computer Programming, Volume 2: Seminumerical
    // Algorithms, 3rd Edition, 1997, Addison-Wesley, p. 122
    //

    float s;

    do {
        r1 = 2 * uniform(state) - 1;
        r2 = 2 * uniform(state) - 1;
        s = r1 * r1 + r2 * r2;
    } while (s >= 1);

    s = sqrtf(-2 * var * logf(s) / s);
    r1 *= s;
    r2 *= s;
}

__device__ void normal(float4& v, float var, state_type& state)
{
    normal(v.x, v.y, var, state);
    normal(v.z, v.w, var, state);
}

__device__ void normal(float2& v, float var, state_type& state)
{
    normal(v.x, v.y, var, state);
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

/**
 * device function wrappers
 */
cuda::function <void (uint48*)>
    rand48_wrapper::leapfrog(gpu::leapfrog);
cuda::function<void (uint48 const*, uint48 const*, uint48*, uint48*, uint)>
    rand48_wrapper::set(gpu::set);
cuda::function<void (uint48 const*, uint48 const*, uint48*, uint48*, ushort3)>
    rand48_wrapper::restore(gpu::restore);
cuda::function<void (ushort3*)>
    rand48_wrapper::save(gpu::save);
cuda::function<void (float*, uint)>
    rand48_wrapper::uniform(gpu::uniform);
cuda::function<void (uint*, uint)>
    rand48_wrapper::get(gpu::get);

/**
 * device constant wrappers
 */
cuda::symbol<uint48>
    rand48_wrapper::a(gpu::a);
cuda::symbol<uint48>
    rand48_wrapper::c(gpu::c);
cuda::symbol<ushort3*>
    rand48_wrapper::state(gpu::g_state);

}} // namespace random::gpu

} // namespace halmd
