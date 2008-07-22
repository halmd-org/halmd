/* Parallelized rand48 random number generator for CUDA
 *
 * Copyright (C) 2007  Peter Colberg
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

#ifndef MDSIM_GPU_RAND48_H
#define MDSIM_GPU_RAND48_H

#include "cutil.h"


namespace rand48
{

/** leapfrogging multiplier */
extern __constant__ uint3 a;
/** leapfrogging addend */
extern __constant__ uint3 c;


/**
 * combined multiply-add-operation for 48-bit integers
 */
template <typename T>
__inline__ __device__ T muladd(uint3 a, T b, uint3 c)
{
    uint d;
    T r;

    //
    // With a C compiler following the ISO/IEC 9899:1999 standard
    // (described in section 6.3 Conversions), multiplying two 16-bit
    // integers results in promotion to 32-bit integers before the
    // multiplication. The CUDA PTX compiler, however, does not perform
    // integer promotion. Therefore, we have to force conversion to
    // 32-bit integers by using a 32-bit integer type for at least one
    // of the operands.
    //
    // Note that in order to enforce such an integer promotion, it does
    // *not* suffice to add a 32-bit register variable and then assign
    // the value of a 16-bit register to it, as the CUDA compiler will
    // rigorously optimize this additional register away. If in doubt,
    // take a look at the resulting PTX assembly code.
    //

    d = a.x * b.x + c.x;
    r.x = (d & 0xFFFF);

    d >>= 16;

    // Although the next line may overflow we only need the top 16 bits
    // in the following stage, so it does not matter.

    d += a.x * b.y + a.y * b.x + c.y;
    r.y = (d & 0xFFFF);

    d >>= 16;

    d += a.x * b.z + a.y * b.y + a.z * b.x + c.z;
    r.z = (d & 0xFFFF);

    return r;
}


/**
 * returns uniform random number between [0.0, 1.0)
 */
__inline__ __device__ float uniform(ushort3& state)
{
    float r = state.z / 65536.f + state.y / 4294967296.f;
    state = muladd(a, state, c);
    return r;
}


/**
 * generate 2 random numbers from Gaussian distribution with given variance
 */
__inline__ __device__ void gaussian(float& r1, float& r2, float const& var, ushort3& state)
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

} // namespace rand48

#endif /* ! MDSIM_GPU_RAND48_H */
