/* 3-dimensional floating-point vector operations for CUDA device functions
 *
 * Copyright © 2008-2009  Peter Colberg
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

#ifndef HALMD_MATH_GPU_VECTOR3D_CUH
#define HALMD_MATH_GPU_VECTOR3D_CUH

#include <cuda_runtime.h>

namespace halmd { namespace cu
{

template <typename T, unsigned int dimension>
struct vector;

template <>
struct vector<float, 3>
{
    enum { static_size = 3 };
    typedef float value_type;

    float x, y, z;

    __device__ inline vector() {}

    __device__ inline vector(float s) : x(s), y(s), z(s) {}

    __device__ inline vector(float x, float y, float z) : x(x), y(y), z(z) {}

    __device__ inline vector(float3 const& v) : x(v.x), y(v.y), z(v.z) {}

    __device__ inline vector(float4 const& v) : x(v.x), y(v.y), z(v.z) {}

    __device__ inline operator float3() const
    {
        return make_float3(x, y, z);
    }

    __device__ inline operator float4() const
    {
        return make_float4(x, y, z, 0);
    }
};

#ifdef __CUDACC__

/**
 * assignment by componentwise vector<float, 3> addition
 */
__device__ inline vector<float, 3>& operator+=(vector<float, 3>& v, vector<float, 3> const& w)
{
    v.x += w.x;
    v.y += w.y;
    v.z += w.z;
    return v;
}

/**
 * assignment by componentwise vector<float, 3> subtraction
 */
__device__ inline vector<float, 3>& operator-=(vector<float, 3>& v, vector<float, 3> const& w)
{
    v.x -= w.x;
    v.y -= w.y;
    v.z -= w.z;
    return v;
}

/**
 * assignment by scalar multiplication
 */
__device__ inline vector<float, 3>& operator*=(vector<float, 3>& v, float s)
{
    v.x *= s;
    v.y *= s;
    v.z *= s;
    return v;
}

/**
 * assignment by scalar division
 */
__device__ inline vector<float, 3>& operator/=(vector<float, 3>& v, float s)
{
    v.x /= s;
    v.y /= s;
    v.z /= s;
    return v;
}

/**
 * componentwise vector<float, 3> addition
 */
__device__ inline vector<float, 3> operator+(vector<float, 3> v, vector<float, 3> const& w)
{
    v += w;
    return v;
}

/**
 * componentwise vector<float, 3> subtraction
 */
__device__ inline vector<float, 3> operator-(vector<float, 3> v, vector<float, 3> const& w)
{
    v -= w;
    return v;
}

/**
 * componentwise change of sign
 */
__device__ inline vector<float, 3> operator-(vector<float, 3> v)
{
    v.x = -v.x;
    v.y = -v.y;
    v.z = -v.z;
    return v;
}

/**
 * scalar product
 */
float __device__ inline operator*(vector<float, 3> const& v, vector<float, 3> const& w)
{
    return v.x * w.x + v.y * w.y + v.z * w.z;
}

/**
 * scalar multiplication
 */
__device__ inline vector<float, 3> operator*(vector<float, 3> v, float s)
{
    v.x *= s;
    v.y *= s;
    v.z *= s;
    return v;
}

/**
 * scalar multiplication
 */
__device__ inline vector<float, 3> operator*(float s, vector<float, 3> v)
{
    v.x *= s;
    v.y *= s;
    v.z *= s;
    return v;
}

/**
 * scalar division
 */
__device__ inline vector<float, 3> operator/(vector<float, 3> v, float s)
{
    v.x /= s;
    v.y /= s;
    v.z /= s;
    return v;
}

// import predefined CUDA math functions into namespace
using ::ceilf;
using ::cosf;
using ::floorf;
using ::fmodf;
using ::remainderf;
using ::rintf;
using ::roundf;
using ::sinf;
using ::sqrtf;
using ::fabsf;
using ::__fdividef;
using ::__float2int_rd;
using ::__float2int_rn;
using ::__float2int_ru;
using ::__float2int_rz;
using ::__float2uint_rd;
using ::__float2uint_rn;
using ::__float2uint_ru;
using ::__float2uint_rz;
using ::__saturatef;

/**
 * componentwise round to nearest integer
 */
__device__ inline vector<float, 3> rintf(vector<float, 3> v)
{
    v.x = rintf(v.x);
    v.y = rintf(v.y);
    v.z = rintf(v.z);
    return v;
}

/**
 * componentwise round to nearest integer, away from zero
 */
__device__ inline vector<float, 3> roundf(vector<float, 3> v)
{
    v.x = roundf(v.x);
    v.y = roundf(v.y);
    v.z = roundf(v.z);
    return v;
}

/**
 * componentwise round to nearest integer not greater than argument
 */
__device__ inline vector<float, 3> floorf(vector<float, 3> v)
{
    v.x = floorf(v.x);
    v.y = floorf(v.y);
    v.z = floorf(v.z);
    return v;
}

/**
 * componentwise round to nearest integer not less argument
 */
__device__ inline vector<float, 3> ceilf(vector<float, 3> v)
{
    v.x = ceilf(v.x);
    v.y = ceilf(v.y);
    v.z = ceilf(v.z);
    return v;
}

/**
 * componentwise square root function
 */
__device__ inline vector<float, 3> sqrtf(vector<float, 3> v)
{
    v.x = sqrtf(v.x);
    v.y = sqrtf(v.y);
    v.z = sqrtf(v.z);
    return v;
}

/**
 * componentwise cosine function
 */
__device__ inline vector<float, 3> cosf(vector<float, 3> v)
{
    v.x = cosf(v.x);
    v.y = cosf(v.y);
    v.z = cosf(v.z);
    return v;
}

/**
 * componentwise sine function
 */
__device__ inline vector<float, 3> sinf(vector<float, 3> v)
{
    v.x = sinf(v.x);
    v.y = sinf(v.y);
    v.z = sinf(v.z);
    return v;
}

/**
 * componentwise absolute value
 */
__device__ inline vector<float, 3> fabsf(vector<float, 3> v)
{
    v.x = fabsf(v.x);
    v.y = fabsf(v.y);
    v.z = fabsf(v.z);
    return v;
}

/**
 * convert floating-point components to integers, rounding to nearest even integer
 */
__device__ inline int3 __float2int_rn(vector<float, 3> const& v)
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
__device__ inline int3 __float2int_rz(vector<float, 3> const& v)
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
__device__ inline int3 __float2int_ru(vector<float, 3> const& v)
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
__device__ inline int3 __float2int_rd(vector<float, 3> const& v)
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
__device__ inline uint3 __float2uint_rn(vector<float, 3> const& v)
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
__device__ inline uint3 __float2uint_rz(vector<float, 3> const& v)
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
__device__ inline uint3 __float2uint_ru(vector<float, 3> const& v)
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
__device__ inline uint3 __float2uint_rd(vector<float, 3> const& v)
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
__device__ inline vector<float, 3> __saturatef(vector<float, 3> v)
{
    v.x = __saturatef(v.x);
    v.y = __saturatef(v.y);
    v.z = __saturatef(v.z);
    return v;
}

/**
 * floating-point remainder function, round towards nearest integer
 */
__device__ inline vector<float, 3> remainderf(vector<float, 3> v, float s)
{
    v.x = remainderf(v.x, s);
    v.y = remainderf(v.y, s);
    v.z = remainderf(v.z, s);
    return v;
}

/**
 * floating-point remainder function, round towards zero
 */
__device__ inline vector<float, 3> fmodf(vector<float, 3> v, float s)
{
    v.x = fmodf(v.x, s);
    v.y = fmodf(v.y, s);
    v.z = fmodf(v.z, s);
    return v;
}

/**
 * fast, accurate floating-point division by s < 2^126
 */
__device__ inline vector<float, 3> __fdividef(vector<float, 3> v, float s)
{
    v.x = __fdividef(v.x, s);
    v.y = __fdividef(v.y, s);
    v.z = __fdividef(v.z, s);
    return v;
}

#endif /* __CUDACC__ */

}} // namespace halmd::cu

#endif /* ! HALMD_MATH_GPU_VECTOR3D_CUH */
