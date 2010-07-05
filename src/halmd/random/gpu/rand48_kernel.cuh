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

#ifndef HALMD_RANDOM_GPU_RAND48_KERNEL_CUH
#define HALMD_RANDOM_GPU_RAND48_KERNEL_CUH

#include <halmd/numeric/gpu/uint48.cuh>

namespace halmd
{
namespace random { namespace gpu
{

using numeric::gpu::uint48;

struct rand48_rng
{
    /** per-thread generator state type */
    typedef ushort3 state_type;

    state_type& operator[](unsigned int thread) const
    {
        return g_state[thread];
    }

    /** leapfrogging multiplier */
    uint48 a;
    /** leapfrogging addend */
    uint48 c;
    /** generator states in global device memory */
    state_type* g_state;
};

/**
 * returns uniform random number in [0.0, 1.0)
 */
inline __device__ float uniform(rand48_rng const& rng, rand48_rng::state_type& state)
{
    uint48 const a = rng.a, c = rng.c;
    float variate = state.z / 65536.f + state.y / 4294967296.f;
    state = muladd(a, state, c);
    return variate;
}

/**
 * returns random integer in [0, 2^32-1]
 */
inline __device__ unsigned int get(rand48_rng const& rng, rand48_rng::state_type& state)
{
    uint48 const a = rng.a, c = rng.c;
    unsigned int variate = (state.z << 16UL) + state.y;
    state = muladd(a, state, c);
    return variate;
}

}} // namespace random::gpu

} // namespace halmd

#endif /* ! HALMD_RANDOM_GPU_RAND48_KERNEL_CUH */
