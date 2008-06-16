/* 2-dimensional floating-point vector operations for CUDA device functions
 *
 * Copyright (C) 2008  Peter Colberg
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

#ifndef MDSIM_GPU_VECTOR2D_H
#define MDSIM_GPU_VECTOR2D_H

#include <cuda/cuda_runtime.h>


/**
 * equality comparison
 */
__device__ bool operator==(float2 const& v, float2 const& w)
{
    return (v.x == w.x && v.y == w.y);
}

/**
 * inequality comparison
 */
__device__ bool operator!=(float2 const& v, float2 const& w)
{
    return (v.x != w.x || v.y != w.y);
}

/**
 * componentwise less than comparison
 */
__device__ bool operator<(float2 const& v, float2 const& w)
{
    return (v.x < w.x && v.y < w.y);
}

/**
 * componentwise greater than comparison
 */
__device__ bool operator>(float2 const& v, float2 const& w)
{
    return (v.x > w.x && v.y > w.y);
}

/**
 * componentwise less than or equal to comparison
 */
__device__ bool operator<=(float2 const& v, float2 const& w)
{
    return (v.x <= w.x && v.y <= w.y);
}

/**
 * componentwise greater than or equal to comparison
 */
__device__ bool operator>=(float2 const& v, float2 const& w)
{
    return (v.x >= w.x && v.y >= w.y);
}

/**
 * assignment by componentwise vector addition
 */
__device__ float2& operator+=(float2& v, float2 const& w)
{
    v.x += w.x;
    v.y += w.y;
    return v;
}

/**
 * assignment by componentwise vector subtraction
 */
__device__ float2& operator-=(float2& v, float2 const& w)
{
    v.x -= w.x;
    v.y -= w.y;
    return v;
}

/**
 * assignment by scalar multiplication
 */
__device__ float2& operator*=(float2& v, float const& s)
{
    v.x *= s;
    v.y *= s;
    return v;
}

/**
 * assignment by scalar division
 */
__device__ float2& operator/=(float2& v, float const& s)
{
    v.x /= s;
    v.y /= s;
    return v;
}

/**
 * componentwise vector addition
 */
__device__ float2 operator+(float2 v, float2 const& w)
{
    v.x += w.x;
    v.y += w.y;
    return v;
}

/**
 * componentwise vector subtraction
 */
__device__ float2 operator-(float2 v, float2 const& w)
{
    v.x -= w.x;
    v.y -= w.y;
    return v;
}

/**
 * scalar product
 */
__device__ float operator*(float2 const& v, float2 const& w)
{
    return v.x * w.x + v.y * w.y;
}

/**
 * scalar multiplication
 */
__device__ float2 operator*(float2 v, float const& s)
{
    v.x *= s;
    v.y *= s;
    return v;
}

/**
 * scalar multiplication
 */
__device__ float2 operator*(float const& s, float2 v)
{
    v.x *= s;
    v.y *= s;
    return v;
}

/**
 * scalar division
 */
__device__ float2 operator/(float2 v, float const& s)
{
    v.x /= s;
    v.y /= s;
    return v;
}

/**
 * componentwise round to nearest integer
 */
__device__ float2 rintf(float2 v)
{
    v.x = rintf(v.x);
    v.y = rintf(v.y);
    return v;
}


/**
 * componentwise round to nearest integer, away from zero
 */
__device__ float2 roundf(float2 v)
{
    v.x = roundf(v.x);
    v.y = roundf(v.y);
    return v;
}


/**
 * componentwise round to nearest integer not greater than argument
 */
__device__ float2 floorf(float2 v)
{
    v.x = floorf(v.x);
    v.y = floorf(v.y);
    return v;
}


/**
 * componentwise round to nearest integer not less argument
 */
__device__ float2 ceilf(float2 v)
{
    v.x = ceilf(v.x);
    v.y = ceilf(v.y);
    return v;
}

/**
 * componentwise square root function
 */
__device__ float2 sqrtf(float2 v)
{
    v.x = sqrtf(v.x);
    v.y = sqrtf(v.y);
    return v;
}


/**
 * convert floating-point components to integers, rounding to nearest even integer
 */
__device__ int2 __float2int_rn(float2 const& v)
{
    int2 w;
    w.x = __float2int_rn(v.x);
    w.y = __float2int_rn(v.y);
    return w;
}

/**
 * convert floating-point components to integers, rounding towards zero
 */
__device__ int2 __float2int_rz(float2 const& v)
{
    int2 w;
    w.x = __float2int_rz(v.x);
    w.y = __float2int_rz(v.y);
    return w;
}

/**
 * convert floating-point components to integers, rounding to positive infinity
 */
__device__ int2 __float2int_ru(float2 const& v)
{
    int2 w;
    w.x = __float2int_ru(v.x);
    w.y = __float2int_ru(v.y);
    return w;
}

/**
 * convert floating-point components to integers, rounding to negative infinity
 */
__device__ int2 __float2int_rd(float2 const& v)
{
    int2 w;
    w.x = __float2int_rd(v.x);
    w.y = __float2int_rd(v.y);
    return w;
}

/**
 * convert floating-point components to unsigned integers, rounding to nearest even integer
 */
__device__ uint2 __float2uint_rn(float2 const& v)
{
    uint2 w;
    w.x = __float2uint_rn(v.x);
    w.y = __float2uint_rn(v.y);
    return w;
}

/**
 * convert floating-point components to unsigned integers, rounding towards zero
 */
__device__ uint2 __float2uint_rz(float2 const& v)
{
    uint2 w;
    w.x = __float2uint_rz(v.x);
    w.y = __float2uint_rz(v.y);
    return w;
}

/**
 * convert floating-point components to unsigned integers, rounding to positive infinity
 */
__device__ uint2 __float2uint_ru(float2 const& v)
{
    uint2 w;
    w.x = __float2uint_ru(v.x);
    w.y = __float2uint_ru(v.y);
    return w;
}

/**
 * convert floating-point components to unsigned integers, rounding to negative infinity
 */
__device__ uint2 __float2uint_rd(float2 const& v)
{
    uint2 w;
    w.x = __float2uint_rd(v.x);
    w.y = __float2uint_rd(v.y);
    return w;
}

/**
 * limit floating-point components to unit interval [0, 1]
 */
__device__ float2 __saturatef(float2 v)
{
    v.x = __saturatef(v.x);
    v.y = __saturatef(v.y);
    return v;
}

/**
 * floating-point remainder function, round towards nearest integer
 */
__device__ float2 remainderf(float2 v, const float s)
{
    v.x = remainderf(v.x, s);
    v.y = remainderf(v.y, s);
    return v;
}

/**
 * floating-point remainder function, round towards zero
 */
__device__ float2 fmodf(float2 v, const float s)
{
    v.x = fmodf(v.x, s);
    v.y = fmodf(v.y, s);
    return v;
}

/**
 * fast, accurate floating-point division by s < 2^126
 */
__device__ float2 __fdividef(float2 v, const float s)
{
    v.x = __fdividef(v.x, s);
    v.y = __fdividef(v.y, s);
    return v;
}

#endif /* ! MDSIM_GPU_VECTOR2D_H */
