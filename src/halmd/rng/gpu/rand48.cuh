/* Parallelized rand48 random number generator for CUDA
 *
 * Copyright © 2007-2009  Peter Colberg
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

#ifndef HALMD_RNG_GPU_RAND48_H
#define HALMD_RNG_GPU_RAND48_H

#include <halmd/rng/gpu/uint48.cuh>
#include <halmd/math/gpu/vector2d.cuh>
#include <halmd/math/gpu/vector3d.cuh>

namespace halmd { namespace cu
{
#ifdef CU_NAMESPACE
namespace CU_NAMESPACE
{
#endif
namespace rand48
{

typedef ushort3 state_type;

/** leapfrogging multiplier */
__constant__ uint48 a;
/** leapfrogging addend */
__constant__ uint48 c;
/** generator state in global device memory */
__constant__ state_type* g_state;

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
__device__ void gaussian(float& r1, float& r2, float var, state_type& state)
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

__device__ void gaussian(float4& v, float var, state_type& state)
{
    gaussian(v.x, v.y, var, state);
    gaussian(v.z, v.w, var, state);
}

__device__ void gaussian(float2& v, float var, state_type& state)
{
    gaussian(v.x, v.y, var, state);
}

}
#ifdef CU_NAMESPACE
}
#endif
}}

#endif /* ! HALMD_RNG_GPU_RAND48_H */
