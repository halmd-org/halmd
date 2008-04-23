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

#endif /* ! MDSIM_GPU_VECTOR2D_H */
