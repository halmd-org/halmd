/* 2-dimensional floating-point vector operations for CUDA device functions
 *
 * Copyright © 2008-2009  Peter Colberg
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

#ifndef LJGPU_MATH_GPU_VECTOR2D_CUH
#define LJGPU_MATH_GPU_VECTOR2D_CUH

#include <cuda_runtime.h>

namespace ljgpu { namespace cu
{

template <typename T, unsigned int dimension>
struct vector;

template <>
struct vector<float, 2>
{
    enum { static_size = 2 };
    typedef float value_type;

    float x, y;

    __device__ inline vector() {}

    __device__ inline vector(float s) : x(s), y(s) {}

    __device__ inline vector(float x, float y) : x(x), y(y) {}

    __device__ inline vector(float2 const& v) : x(v.x), y(v.y) {}

    __device__ inline vector(float3 const& v) : x(v.x), y(v.y) {}

    __device__ inline vector(float4 const& v) : x(v.x), y(v.y) {}

    __device__ inline operator float2() const
    {
        return make_float2(x, y);
    }

    __device__ inline operator float3() const
    {
        return make_float3(x, y, 0);
    }

    __device__ inline operator float4() const
    {
        return make_float4(x, y, 0, 0);
    }
};

#ifdef __CUDACC__

/**
 * assignment by componentwise vector<float, 2> addition
 */
__device__ inline vector<float, 2>& operator+=(vector<float, 2>& v, vector<float, 2> const& w)
{
    v.x += w.x;
    v.y += w.y;
    return v;
}

/**
 * assignment by componentwise vector<float, 2> subtraction
 */
__device__ inline vector<float, 2>& operator-=(vector<float, 2>& v, vector<float, 2> const& w)
{
    v.x -= w.x;
    v.y -= w.y;
    return v;
}

/**
 * assignment by scalar multiplication
 */
__device__ inline vector<float, 2>& operator*=(vector<float, 2>& v, float s)
{
    v.x *= s;
    v.y *= s;
    return v;
}

/**
 * assignment by scalar division
 */
__device__ inline vector<float, 2>& operator/=(vector<float, 2>& v, float s)
{
    v.x /= s;
    v.y /= s;
    return v;
}

/**
 * componentwise vector<float, 2> addition
 */
__device__ inline vector<float, 2> operator+(vector<float, 2> v, vector<float, 2> const& w)
{
    v += w;
    return v;
}

/**
 * componentwise vector<float, 2> subtraction
 */
__device__ inline vector<float, 2> operator-(vector<float, 2> v, vector<float, 2> const& w)
{
    v -= w;
    return v;
}

/**
 * componentwise change of sign
 */
__device__ inline vector<float, 2> operator-(vector<float, 2> v)
{
    v.x = -v.x;
    v.y = -v.y;
    return v;
}

/**
 * scalar product
 */
float __device__ inline operator*(vector<float, 2> const& v, vector<float, 2> const& w)
{
    return v.x * w.x + v.y * w.y;
}

/**
 * scalar multiplication
 */
__device__ inline vector<float, 2> operator*(vector<float, 2> v, float s)
{
    v.x *= s;
    v.y *= s;
    return v;
}

/**
 * scalar multiplication
 */
__device__ inline vector<float, 2> operator*(float s, vector<float, 2> v)
{
    v.x *= s;
    v.y *= s;
    return v;
}

/**
 * scalar division
 */
__device__ inline vector<float, 2> operator/(vector<float, 2> v, float s)
{
    v.x /= s;
    v.y /= s;
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
__device__ inline vector<float, 2> rintf(vector<float, 2> v)
{
    v.x = rintf(v.x);
    v.y = rintf(v.y);
    return v;
}

/**
 * componentwise round to nearest integer, away from zero
 */
__device__ inline vector<float, 2> roundf(vector<float, 2> v)
{
    v.x = roundf(v.x);
    v.y = roundf(v.y);
    return v;
}

/**
 * componentwise round to nearest integer not greater than argument
 */
__device__ inline vector<float, 2> floorf(vector<float, 2> v)
{
    v.x = floorf(v.x);
    v.y = floorf(v.y);
    return v;
}

/**
 * componentwise round to nearest integer not less argument
 */
__device__ inline vector<float, 2> ceilf(vector<float, 2> v)
{
    v.x = ceilf(v.x);
    v.y = ceilf(v.y);
    return v;
}

/**
 * componentwise square root function
 */
__device__ inline vector<float, 2> sqrtf(vector<float, 2> v)
{
    v.x = sqrtf(v.x);
    v.y = sqrtf(v.y);
    return v;
}

/**
 * componentwise cosine function
 */
__device__ inline vector<float, 2> cosf(vector<float, 2> v)
{
    v.x = cosf(v.x);
    v.y = cosf(v.y);
    return v;
}

/**
 * componentwise sine function
 */
__device__ inline vector<float, 2> sinf(vector<float, 2> v)
{
    v.x = sinf(v.x);
    v.y = sinf(v.y);
    return v;
}

/**
 * componentwise absolute value
 */
__device__ inline vector<float, 2> fabsf(vector<float, 2> v)
{
    v.x = fabsf(v.x);
    v.y = fabsf(v.y);
    return v;
}

/**
 * convert floating-point components to integers, rounding to nearest even integer
 */
__device__ inline int2 __float2int_rn(vector<float, 2> const& v)
{
    int2 w;
    w.x = __float2int_rn(v.x);
    w.y = __float2int_rn(v.y);
    return w;
}

/**
 * convert floating-point components to integers, rounding towards zero
 */
__device__ inline int2 __float2int_rz(vector<float, 2> const& v)
{
    int2 w;
    w.x = __float2int_rz(v.x);
    w.y = __float2int_rz(v.y);
    return w;
}

/**
 * convert floating-point components to integers, rounding to positive infinity
 */
__device__ inline int2 __float2int_ru(vector<float, 2> const& v)
{
    int2 w;
    w.x = __float2int_ru(v.x);
    w.y = __float2int_ru(v.y);
    return w;
}

/**
 * convert floating-point components to integers, rounding to negative infinity
 */
__device__ inline int2 __float2int_rd(vector<float, 2> const& v)
{
    int2 w;
    w.x = __float2int_rd(v.x);
    w.y = __float2int_rd(v.y);
    return w;
}

/**
 * convert floating-point components to unsigned integers, rounding to nearest even integer
 */
__device__ inline uint2 __float2uint_rn(vector<float, 2> const& v)
{
    uint2 w;
    w.x = __float2uint_rn(v.x);
    w.y = __float2uint_rn(v.y);
    return w;
}

/**
 * convert floating-point components to unsigned integers, rounding towards zero
 */
__device__ inline uint2 __float2uint_rz(vector<float, 2> const& v)
{
    uint2 w;
    w.x = __float2uint_rz(v.x);
    w.y = __float2uint_rz(v.y);
    return w;
}

/**
 * convert floating-point components to unsigned integers, rounding to positive infinity
 */
__device__ inline uint2 __float2uint_ru(vector<float, 2> const& v)
{
    uint2 w;
    w.x = __float2uint_ru(v.x);
    w.y = __float2uint_ru(v.y);
    return w;
}

/**
 * convert floating-point components to unsigned integers, rounding to negative infinity
 */
__device__ inline uint2 __float2uint_rd(vector<float, 2> const& v)
{
    uint2 w;
    w.x = __float2uint_rd(v.x);
    w.y = __float2uint_rd(v.y);
    return w;
}

/**
 * limit floating-point components to unit interval [0, 1]
 */
__device__ inline vector<float, 2> __saturatef(vector<float, 2> v)
{
    v.x = __saturatef(v.x);
    v.y = __saturatef(v.y);
    return v;
}

/**
 * floating-point remainder function, round towards nearest integer
 */
__device__ inline vector<float, 2> remainderf(vector<float, 2> v, float s)
{
    v.x = remainderf(v.x, s);
    v.y = remainderf(v.y, s);
    return v;
}

/**
 * floating-point remainder function, round towards zero
 */
__device__ inline vector<float, 2> fmodf(vector<float, 2> v, float s)
{
    v.x = fmodf(v.x, s);
    v.y = fmodf(v.y, s);
    return v;
}

/**
 * fast, accurate floating-point division by s < 2^126
 */
__device__ inline vector<float, 2> __fdividef(vector<float, 2> v, float s)
{
    v.x = __fdividef(v.x, s);
    v.y = __fdividef(v.y, s);
    return v;
}

#endif /* __CUDACC__ */

}} // namespace ljgpu::cu

#endif /* ! LJGPU_MATH_GPU_VECTOR2D_CUH */
