/* 4-dimensional floating-point vector operations for CUDA device functions
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

#ifndef HALMD_MATH_GPU_VECTOR4D_CUH
#define HALMD_MATH_GPU_VECTOR4D_CUH

#include <cuda_runtime.h>

namespace halmd { namespace cu
{

// To prevent shadowing the math functions of the :: namespace,
// we define the vector class in a child namespace and later
// import the vector class into the parent namespace. Through the
// magic of Koenig lookup, the vector math functions defined in
// the child namespace will still be found by the compiler.
namespace detail
{

template <typename T, unsigned int dimension>
struct vector;

template <>
struct vector<float, 4>
{
    enum { static_size = 4 };
    typedef float value_type;

    float x, y, z, w;

    __device__ inline vector() {}

    __device__ inline vector(float s) : x(s), y(s), z(s), w(s) {}

    __device__ inline vector(float x, float y, float z, float w) : x(x), y(y), z(z), w(w) {}

    __device__ inline vector(float4 const& v) : x(v.x), y(v.y), z(v.z), w(v.w) {}

    __device__ inline operator float4() const
    {
        return make_float4(x, y, z, w);
    }
};

#ifdef __CUDACC__

/**
 * assignment by componentwise vector<float, 4> addition
 */
__device__ inline vector<float, 4>& operator+=(vector<float, 4>& v, vector<float, 4> const& w)
{
    v.x += w.x;
    v.y += w.y;
    v.z += w.z;
    v.w += w.w;
    return v;
}

/**
 * assignment by componentwise vector<float, 4> subtraction
 */
__device__ inline vector<float, 4>& operator-=(vector<float, 4>& v, vector<float, 4> const& w)
{
    v.x -= w.x;
    v.y -= w.y;
    v.z -= w.z;
    v.w -= w.w;
    return v;
}

/**
 * assignment by scalar multiplication
 */
__device__ inline vector<float, 4>& operator*=(vector<float, 4>& v, float s)
{
    v.x *= s;
    v.y *= s;
    v.z *= s;
    v.w *= s;
    return v;
}

/**
 * assignment by scalar division
 */
__device__ inline vector<float, 4>& operator/=(vector<float, 4>& v, float s)
{
    v.x /= s;
    v.y /= s;
    v.z /= s;
    v.w /= s;
    return v;
}

/**
 * componentwise vector<float, 4> addition
 */
__device__ inline vector<float, 4> operator+(vector<float, 4> v, vector<float, 4> const& w)
{
    v += w;
    return v;
}

/**
 * componentwise vector<float, 4> subtraction
 */
__device__ inline vector<float, 4> operator-(vector<float, 4> v, vector<float, 4> const& w)
{
    v -= w;
    return v;
}

/**
 * componentwise change of sign
 */
__device__ inline vector<float, 4> operator-(vector<float, 4> v)
{
    v.x = -v.x;
    v.y = -v.y;
    v.z = -v.z;
    v.w = -v.w;
    return v;
}

/**
 * scalar product
 */
float __device__ inline operator*(vector<float, 4> const& v, vector<float, 4> const& w)
{
    return v.x * w.x + v.y * w.y + v.z * w.z + v.w * w.w;
}

/**
 * scalar multiplication
 */
__device__ inline vector<float, 4> operator*(vector<float, 4> v, float s)
{
    v.x *= s;
    v.y *= s;
    v.z *= s;
    v.w *= s;
    return v;
}

/**
 * scalar multiplication
 */
__device__ inline vector<float, 4> operator*(float s, vector<float, 4> v)
{
    v.x *= s;
    v.y *= s;
    v.z *= s;
    v.w *= s;
    return v;
}

/**
 * scalar division
 */
__device__ inline vector<float, 4> operator/(vector<float, 4> v, float s)
{
    v.x /= s;
    v.y /= s;
    v.z /= s;
    v.w /= s;
    return v;
}

/**
 * componentwise round to nearest integer
 */
__device__ inline vector<float, 4> rintf(vector<float, 4> v)
{
    v.x = ::rintf(v.x);
    v.y = ::rintf(v.y);
    v.z = ::rintf(v.z);
    v.w = ::rintf(v.w);
    return v;
}

/**
 * componentwise round to nearest integer, away from zero
 */
__device__ inline vector<float, 4> roundf(vector<float, 4> v)
{
    v.x = ::roundf(v.x);
    v.y = ::roundf(v.y);
    v.z = ::roundf(v.z);
    v.w = ::roundf(v.w);
    return v;
}

/**
 * componentwise round to nearest integer not greater than argument
 */
__device__ inline vector<float, 4> floorf(vector<float, 4> v)
{
    v.x = ::floorf(v.x);
    v.y = ::floorf(v.y);
    v.z = ::floorf(v.z);
    v.w = ::floorf(v.w);
    return v;
}

/**
 * componentwise round to nearest integer not less argument
 */
__device__ inline vector<float, 4> ceilf(vector<float, 4> v)
{
    v.x = ::ceilf(v.x);
    v.y = ::ceilf(v.y);
    v.z = ::ceilf(v.z);
    v.w = ::ceilf(v.w);
    return v;
}

/**
 * componentwise square root function
 */
__device__ inline vector<float, 4> sqrtf(vector<float, 4> v)
{
    v.x = ::sqrtf(v.x);
    v.y = ::sqrtf(v.y);
    v.z = ::sqrtf(v.z);
    v.w = ::sqrtf(v.w);
    return v;
}

/**
 * componentwise cosine function
 */
__device__ inline vector<float, 4> cosf(vector<float, 4> v)
{
    v.x = ::cosf(v.x);
    v.y = ::cosf(v.y);
    v.z = ::cosf(v.z);
    v.w = ::cosf(v.w);
    return v;
}

/**
 * componentwise sine function
 */
__device__ inline vector<float, 4> sinf(vector<float, 4> v)
{
    v.x = ::sinf(v.x);
    v.y = ::sinf(v.y);
    v.z = ::sinf(v.z);
    v.w = ::sinf(v.w);
    return v;
}

/**
 * componentwise absolute value
 */
__device__ inline vector<float, 4> fabsf(vector<float, 4> v)
{
    v.x = ::fabsf(v.x);
    v.y = ::fabsf(v.y);
    v.z = ::fabsf(v.z);
    v.w = ::fabsf(v.w);
    return v;
}

/**
 * convert floating-point components to integers, rounding to nearest even integer
 */
__device__ inline int4 __float2int_rn(vector<float, 4> const& v)
{
    int4 w;
    w.x = ::__float2int_rn(v.x);
    w.y = ::__float2int_rn(v.y);
    w.z = ::__float2int_rn(v.z);
    w.w = ::__float2int_rn(v.w);
    return w;
}

/**
 * convert floating-point components to integers, rounding towards zero
 */
__device__ inline int4 __float2int_rz(vector<float, 4> const& v)
{
    int4 w;
    w.x = ::__float2int_rz(v.x);
    w.y = ::__float2int_rz(v.y);
    w.z = ::__float2int_rz(v.z);
    w.w = ::__float2int_rz(v.w);
    return w;
}

/**
 * convert floating-point components to integers, rounding to positive infinity
 */
__device__ inline int4 __float2int_ru(vector<float, 4> const& v)
{
    int4 w;
    w.x = ::__float2int_ru(v.x);
    w.y = ::__float2int_ru(v.y);
    w.z = ::__float2int_ru(v.z);
    w.w = ::__float2int_ru(v.w);
    return w;
}

/**
 * convert floating-point components to integers, rounding to negative infinity
 */
__device__ inline int4 __float2int_rd(vector<float, 4> const& v)
{
    int4 w;
    w.x = ::__float2int_rd(v.x);
    w.y = ::__float2int_rd(v.y);
    w.z = ::__float2int_rd(v.z);
    w.w = ::__float2int_rd(v.w);
    return w;
}

/**
 * convert floating-point components to unsigned integers, rounding to nearest even integer
 */
__device__ inline uint4 __float2uint_rn(vector<float, 4> const& v)
{
    uint4 w;
    w.x = ::__float2uint_rn(v.x);
    w.y = ::__float2uint_rn(v.y);
    w.z = ::__float2uint_rn(v.z);
    w.w = ::__float2uint_rn(v.w);
    return w;
}

/**
 * convert floating-point components to unsigned integers, rounding towards zero
 */
__device__ inline uint4 __float2uint_rz(vector<float, 4> const& v)
{
    uint4 w;
    w.x = ::__float2uint_rz(v.x);
    w.y = ::__float2uint_rz(v.y);
    w.z = ::__float2uint_rz(v.z);
    w.w = ::__float2uint_rz(v.w);
    return w;
}

/**
 * convert floating-point components to unsigned integers, rounding to positive infinity
 */
__device__ inline uint4 __float2uint_ru(vector<float, 4> const& v)
{
    uint4 w;
    w.x = ::__float2uint_ru(v.x);
    w.y = ::__float2uint_ru(v.y);
    w.z = ::__float2uint_ru(v.z);
    w.w = ::__float2uint_ru(v.w);
    return w;
}

/**
 * convert floating-point components to unsigned integers, rounding to negative infinity
 */
__device__ inline uint4 __float2uint_rd(vector<float, 4> const& v)
{
    uint4 w;
    w.x = ::__float2uint_rd(v.x);
    w.y = ::__float2uint_rd(v.y);
    w.z = ::__float2uint_rd(v.z);
    w.w = ::__float2uint_rd(v.w);
    return w;
}

/**
 * limit floating-point components to unit interval [0, 1]
 */
__device__ inline vector<float, 4> __saturatef(vector<float, 4> v)
{
    v.x = ::__saturatef(v.x);
    v.y = ::__saturatef(v.y);
    v.z = ::__saturatef(v.z);
    v.w = ::__saturatef(v.w);
    return v;
}

/**
 * floating-point remainder function, round towards nearest integer
 */
__device__ inline vector<float, 4> remainderf(vector<float, 4> v, float s)
{
    v.x = ::remainderf(v.x, s);
    v.y = ::remainderf(v.y, s);
    v.z = ::remainderf(v.z, s);
    v.w = ::remainderf(v.w, s);
    return v;
}

/**
 * floating-point remainder function, round towards zero
 */
__device__ inline vector<float, 4> fmodf(vector<float, 4> v, float s)
{
    v.x = ::fmodf(v.x, s);
    v.y = ::fmodf(v.y, s);
    v.z = ::fmodf(v.z, s);
    v.w = ::fmodf(v.w, s);
    return v;
}

/**
 * fast, accurate floating-point division by s < 2^126
 */
__device__ inline vector<float, 4> __fdividef(vector<float, 4> v, float s)
{
    v.x = ::__fdividef(v.x, s);
    v.y = ::__fdividef(v.y, s);
    v.z = ::__fdividef(v.z, s);
    v.w = ::__fdividef(v.w, s);
    return v;
}

#endif /* __CUDACC__ */

} // namespace detail

// import vector class into parent namespace
using detail::vector;

}} // namespace halmd::cu

#endif /* ! HALMD_MATH_GPU_VECTOR4D_CUH */
