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

#ifndef HALMD_NUMERIC_MP_UINT48_HPP
#define HALMD_NUMERIC_MP_UINT48_HPP

#include <halmd/config.hpp>
#include <halmd/numeric/zero.hpp>

//
// The multiply-add operation is based on the rand48 generator of the
// GNU Scientific Library. The file rng/rand48.c was written by James
// Theiler and Brian Gough. The code snippet in muladd() is merely
// 8 lines long and can thus be used freely.
//

namespace halmd {

/**
 * We define an 48 bit integer type, along with multiply-add and
 * add operators. It is used in the rand48 random number generator
 * for CUDA, which lacks native 64 bit integers.
 *
 * Note that we use three 32 bit integers instead of 16 bit
 * integers, which is reasoned in the multiply-add function,
 * and that we define a custom uint48 type instead of using
 * the predefined uint3, as we overload the add operators.
 */
struct uint48
{
    unsigned int x, y, z;
};

inline HALMD_GPU_ENABLED uint48 make_uint48(unsigned int x, unsigned int y, unsigned int z)
{
    uint48 u; u.x = x; u.y = y; u.z = z; return u;
}

/**
 * Helper to initialize 48-bit integer to zero.
 */
template <>
struct zero<uint48>
{
    HALMD_GPU_ENABLED operator uint48() const
    {
        return make_uint48(0, 0, 0);
    }
};

/**
 * combined multiply-add operation for 48 bit integers
 */
template <typename T>
HALMD_GPU_ENABLED T muladd(uint48 const& a, T const& b, uint48 const& c)
{
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

    unsigned int d = a.x * b.x + c.x;
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
 * add-operation for 48 bit integers
 */
inline HALMD_GPU_ENABLED uint48& operator+=(uint48& a, uint48 const& b)
{
    unsigned int d = a.x + b.x;
    a.x = (d & 0xFFFF);
    d >>= 16;
    d += a.y + b.y;
    a.y = (d & 0xFFFF);
    d >>= 16;
    d += a.z + b.z;
    a.z = (d & 0xFFFF);
    return a;
}

/**
 * add-operation for 48 bit integers
 */
inline HALMD_GPU_ENABLED uint48 operator+(uint48 const& a, uint48 b)
{
    b += a;
    return b;
}

} // namespace halmd

#endif /* ! HALMD_NUMERIC_MP_UINT48_HPP */
