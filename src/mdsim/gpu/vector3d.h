/* 3-dimensional floating-point vector operations for CUDA device functions
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

#ifndef MDSIM_GPU_VECTOR3D_H
#define MDSIM_GPU_VECTOR3D_H

#include <cuda/cuda_runtime.h>


/**
 * equality comparison
 */
__device__ bool operator==(float3 const& v, float3 const& w)
{
    return (v.x == w.x && v.y == w.y && v.z == w.z);
}

/**
 * inequality comparison
 */
__device__ bool operator!=(float3 const& v, float3 const& w)
{
    return (v.x != w.x || v.y != w.y || v.z != w.z);
}

/**
 * componentwise less than comparison
 */
__device__ bool operator<(float3 const& v, float3 const& w)
{
    return (v.x < w.x && v.y < w.y && v.z < w.z);
}

/**
 * componentwise greater than comparison
 */
__device__ bool operator>(float3 const& v, float3 const& w)
{
    return (v.x > w.x && v.y > w.y && v.z > w.z);
}

/**
 * componentwise less than or equal to comparison
 */
__device__ bool operator<=(float3 const& v, float3 const& w)
{
    return (v.x <= w.x && v.y <= w.y && v.z <= w.z);
}

/**
 * componentwise greater than or equal to comparison
 */
__device__ bool operator>=(float3 const& v, float3 const& w)
{
    return (v.x >= w.x && v.y >= w.y && v.z >= w.z);
}

/**
 * assignment by componentwise vector addition
 */
__device__ float3& operator+=(float3& v, float3 const& w)
{
    v.x += w.x;
    v.y += w.y;
    v.z += w.z;
    return v;
}

/**
 * assignment by componentwise vector subtraction
 */
__device__ float3& operator-=(float3& v, float3 const& w)
{
    v.x -= w.x;
    v.y -= w.y;
    v.z -= w.z;
    return v;
}

/**
 * assignment by scalar multiplication
 */
__device__ float3& operator*=(float3& v, float const& s)
{
    v.x *= s;
    v.y *= s;
    v.z *= s;
    return v;
}

/**
 * assignment by scalar division
 */
__device__ float3& operator/=(float3& v, float const& s)
{
    v.x /= s;
    v.y /= s;
    v.z /= s;
    return v;
}

/**
 * componentwise vector addition
 */
__device__ float3 operator+(float3 v, float3 const& w)
{
    v.x += w.x;
    v.y += w.y;
    v.z += w.z;
    return v;
}

/**
 * componentwise vector subtraction
 */
__device__ float3 operator-(float3 v, float3 const& w)
{
    v.x -= w.x;
    v.y -= w.y;
    v.z -= w.z;
    return v;
}

/**
 * scalar product
 */
__device__ float operator*(float3 const& v, float3 const& w)
{
    return v.x * w.x + v.y * w.y + v.z * w.z;
}

/**
 * scalar multiplication
 */
__device__ float3 operator*(float3 v, float const& s)
{
    v.x *= s;
    v.y *= s;
    v.z *= s;
    return v;
}

/**
 * scalar multiplication
 */
__device__ float3 operator*(float const& s, float3 v)
{
    v.x *= s;
    v.y *= s;
    v.z *= s;
    return v;
}

/**
 * scalar division
 */
__device__ float3 operator/(float3 v, float const& s)
{
    v.x /= s;
    v.y /= s;
    v.z /= s;
    return v;
}

/**
 * componentwise round to nearest integer
 */
__device__ float3 rintf(float3 v)
{
    v.x = rintf(v.x);
    v.y = rintf(v.y);
    v.z = rintf(v.z);
    return v;
}


/**
 * componentwise round to nearest integer, away from zero
 */
__device__ float3 roundf(float3 v)
{
    v.x = roundf(v.x);
    v.y = roundf(v.y);
    v.z = roundf(v.z);
    return v;
}


/**
 * componentwise round to nearest integer not greater than argument
 */
__device__ float3 floorf(float3 v)
{
    v.x = floorf(v.x);
    v.y = floorf(v.y);
    v.z = floorf(v.z);
    return v;
}


/**
 * componentwise round to nearest integer not less argument
 */
__device__ float3 ceilf(float3 v)
{
    v.x = ceilf(v.x);
    v.y = ceilf(v.y);
    v.z = ceilf(v.z);
    return v;
}

/**
 * componentwise square root function
 */
__device__ float3 sqrtf(float3 v)
{
    v.x = sqrtf(v.x);
    v.y = sqrtf(v.y);
    v.z = sqrtf(v.z);
    return v;
}


/**
 * convert floating-point components to integers, rounding to nearest even integer
 */
__device__ int3 __float2int_rn(float3 const& v)
{
    int3 w;
    w.x = __float2int_rn(v.x);
    w.y = __float2int_rn(v.y);
    w.z = __float2int_rn(v.z);
    return w;
}

/**
 * convert floating-point components to integers, rounding towards zero
 */
__device__ int3 __float2int_rz(float3 const& v)
{
    int3 w;
    w.x = __float2int_rz(v.x);
    w.y = __float2int_rz(v.y);
    w.z = __float2int_rz(v.z);
    return w;
}

/**
 * convert floating-point components to integers, rounding to positive infinity
 */
__device__ int3 __float2int_ru(float3 const& v)
{
    int3 w;
    w.x = __float2int_ru(v.x);
    w.y = __float2int_ru(v.y);
    w.z = __float2int_ru(v.z);
    return w;
}

/**
 * convert floating-point components to integers, rounding to negative infinity
 */
__device__ int3 __float2int_rd(float3 const& v)
{
    int3 w;
    w.x = __float2int_rd(v.x);
    w.y = __float2int_rd(v.y);
    w.z = __float2int_rd(v.z);
    return w;
}

/**
 * convert floating-point components to unsigned integers, rounding to nearest even integer
 */
__device__ uint3 __float2uint_rn(float3 const& v)
{
    uint3 w;
    w.x = __float2uint_rn(v.x);
    w.y = __float2uint_rn(v.y);
    w.z = __float2uint_rn(v.z);
    return w;
}

/**
 * convert floating-point components to unsigned integers, rounding towards zero
 */
__device__ uint3 __float2uint_rz(float3 const& v)
{
    uint3 w;
    w.x = __float2uint_rz(v.x);
    w.y = __float2uint_rz(v.y);
    w.z = __float2uint_rz(v.z);
    return w;
}

/**
 * convert floating-point components to unsigned integers, rounding to positive infinity
 */
__device__ uint3 __float2uint_ru(float3 const& v)
{
    uint3 w;
    w.x = __float2uint_ru(v.x);
    w.y = __float2uint_ru(v.y);
    w.z = __float2uint_ru(v.z);
    return w;
}

/**
 * convert floating-point components to unsigned integers, rounding to negative infinity
 */
__device__ uint3 __float2uint_rd(float3 const& v)
{
    uint3 w;
    w.x = __float2uint_rd(v.x);
    w.y = __float2uint_rd(v.y);
    w.z = __float2uint_rd(v.z);
    return w;
}

/**
 * limit floating-point components to unit interval [0, 1]
 */
__device__ float3 __saturatef(float3 v)
{
    v.x = __saturatef(v.x);
    v.y = __saturatef(v.y);
    v.z = __saturatef(v.z);
    return v;
}

/**
 * floating-point remainder function, round towards nearest integer
 */
__device__ float3 remainderf(float3 v, const float s)
{
    v.x = remainderf(v.x, s);
    v.y = remainderf(v.y, s);
    v.z = remainderf(v.z, s);
    return v;
}

/**
 * floating-point remainder function, round towards zero
 */
__device__ float3 fmodf(float3 v, const float s)
{
    v.x = fmodf(v.x, s);
    v.y = fmodf(v.y, s);
    v.z = fmodf(v.z, s);
    return v;
}

/**
 * fast, accurate floating-point division by s < 2^126
 */
__device__ float3 __fdividef(float3 v, const float s)
{
    v.x = __fdividef(v.x, s);
    v.y = __fdividef(v.y, s);
    v.z = __fdividef(v.z, s);
    return v;
}

#endif /* ! MDSIM_GPU_VECTOR3D_H */
