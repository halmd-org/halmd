/* double-single arithmetic vector operations for CUDA device functions
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

#ifndef LJGPU_MATH_GPU_DSVECTOR_CUH
#define LJGPU_MATH_GPU_DSVECTOR_CUH

#include <cuda_runtime.h>
#include <ljgpu/math/gpu/dsfun.cuh>
#include <ljgpu/math/gpu/vector2d.cuh>
#include <ljgpu/math/gpu/vector3d.cuh>
#include <ljgpu/math/gpu/vector4d.cuh>

namespace ljgpu { namespace cu
{

/**
 * Two-dimensional double-single floating point vector
 */
template <>
struct vector<dfloat, 2>
{
    enum { static_size = 2 };
    typedef dfloat value_type;

    // Use float2 as high- and low-word to work around an NVCC compiler bug:
    // A __shared__ array may only be declared as a struct type if the
    // struct's members are POD or CUDA types, e.g. we cannot use a custom
    // two- or three-dimensional vector of dfloats.
    float2 __hi, __lo;

    __device__ inline vector() {}

    __device__ inline vector(dfloat s) : __hi(make_float2(s.__hi, s.__hi)), __lo(make_float2(s.__lo, s.__lo)) {}

    __device__ inline vector(float s) : __hi(make_float2(s, s)), __lo(make_float2(0, 0)) {}

    __device__ inline vector(vector<float, 2> v, vector<float, 2> w) : __hi(v), __lo(w) {}

    __device__ inline vector(vector<float, 2> v) : __hi(v), __lo(make_float2(0, 0)) {}

    __device__ inline operator vector<float, 2>() const { return __hi; }
};

/**
 * Three-dimensional double-single floating point vector
 */
template <>
struct vector<dfloat, 3>
{
    enum { static_size = 3 };
    typedef dfloat value_type;

    // Use float3 as high- and low-word to work around an NVCC compiler bug:
    // A __shared__ array may only be declared as a struct type if the
    // struct's members are POD or CUDA types, e.g. we cannot use a custom
    // two- or three-dimensional vector of dfloats.
    float3 __hi, __lo;

    __device__ inline vector() {}

    __device__ inline vector(dfloat s) : __hi(make_float3(s.__hi, s.__hi, s.__hi)), __lo(make_float3(s.__lo, s.__lo, s.__lo)) {}

    __device__ inline vector(float s) : __hi(make_float3(s, s, s)), __lo(make_float3(0, 0, 0)) {}

    __device__ inline vector(vector<float, 3> v, vector<float, 3> w) : __hi(v), __lo(w) {}

    __device__ inline vector(vector<float, 3> v) : __hi(v), __lo(make_float3(0, 0, 0)) {}

    __device__ inline operator vector<float, 3>() const { return __hi; }
};

/**
 * Four-dimensional double-single floating point vector
 */
template <>
struct vector<dfloat, 4>
{
    enum { static_size = 4 };
    typedef dfloat value_type;

    // Use float4 as high- and low-word to work around an NVCC compiler bug:
    // A __shared__ array may only be declared as a struct type if the
    // struct's members are POD or CUDA types, e.g. we cannot use a custom
    // two- or three-dimensional vector of dfloats.
    float4 __hi, __lo;

    __device__ inline vector() {}

    __device__ inline vector(dfloat s) : __hi(make_float4(s.__hi, s.__hi, s.__hi, s.__hi)), __lo(make_float4(s.__lo, s.__lo, s.__lo, s.__lo)) {}

    __device__ inline vector(float s) : __hi(make_float4(s, s, s, s)), __lo(make_float4(0, 0, 0, 0)) {}

    __device__ inline vector(vector<float, 4> v, vector<float, 4> w) : __hi(v), __lo(w) {}

    __device__ inline vector(vector<float, 4> v) : __hi(v), __lo(make_float4(0, 0, 0, 0)) {}

    __device__ inline operator vector<float, 4>() const { return __hi; }
};

/**
 * returns high-word floating point vector
 */
template <unsigned int dimension>
__device__ inline vector<float, dimension> dfloat2hi(vector<dfloat, dimension> const& v)
{
    return v.__hi;
}

/**
 * returns low-word floating point vector
 */
template <unsigned int dimension>
__device__ inline vector<float, dimension> dfloat2lo(vector<dfloat, dimension> const& v)
{
    return v.__lo;
}

/**
 * assignment by componentwise vector addition
 */
__device__ inline vector<dfloat, 2>& operator+=(vector<dfloat, 2>& v, vector<dfloat, 2> const& w)
{
    __dsadd(v.__hi.x, v.__lo.x, v.__hi.x, v.__lo.x, w.__hi.x, w.__lo.x);
    __dsadd(v.__hi.y, v.__lo.y, v.__hi.y, v.__lo.y, w.__hi.y, w.__lo.y);
    return v;
}

__device__ inline vector<dfloat, 3>& operator+=(vector<dfloat, 3>& v, vector<dfloat, 3> const& w)
{
    __dsadd(v.__hi.x, v.__lo.x, v.__hi.x, v.__lo.x, w.__hi.x, w.__lo.x);
    __dsadd(v.__hi.y, v.__lo.y, v.__hi.y, v.__lo.y, w.__hi.y, w.__lo.y);
    __dsadd(v.__hi.z, v.__lo.z, v.__hi.z, v.__lo.z, w.__hi.z, w.__lo.z);
    return v;
}

__device__ inline vector<dfloat, 4>& operator+=(vector<dfloat, 4>& v, vector<dfloat, 4> const& w)
{
    __dsadd(v.__hi.x, v.__lo.x, v.__hi.x, v.__lo.x, w.__hi.x, w.__lo.x);
    __dsadd(v.__hi.y, v.__lo.y, v.__hi.y, v.__lo.y, w.__hi.y, w.__lo.y);
    __dsadd(v.__hi.z, v.__lo.z, v.__hi.z, v.__lo.z, w.__hi.z, w.__lo.z);
    __dsadd(v.__hi.w, v.__lo.w, v.__hi.w, v.__lo.w, w.__hi.w, w.__lo.w);
    return v;
}

/**
 * assignment by componentwise vector subtraction
 */
__device__ inline vector<dfloat, 2>& operator-=(vector<dfloat, 2>& v, vector<dfloat, 2> const& w)
{
    __dssub(v.__hi.x, v.__lo.x, v.__hi.x, v.__lo.x, w.__hi.x, w.__lo.x);
    __dssub(v.__hi.y, v.__lo.y, v.__hi.y, v.__lo.y, w.__hi.y, w.__lo.y);
    return v;
}

__device__ inline vector<dfloat, 3>& operator-=(vector<dfloat, 3>& v, vector<dfloat, 3> const& w)
{
    __dssub(v.__hi.x, v.__lo.x, v.__hi.x, v.__lo.x, w.__hi.x, w.__lo.x);
    __dssub(v.__hi.y, v.__lo.y, v.__hi.y, v.__lo.y, w.__hi.y, w.__lo.y);
    __dssub(v.__hi.z, v.__lo.z, v.__hi.z, v.__lo.z, w.__hi.z, w.__lo.z);
    return v;
}

__device__ inline vector<dfloat, 4>& operator-=(vector<dfloat, 4>& v, vector<dfloat, 4> const& w)
{
    __dssub(v.__hi.x, v.__lo.x, v.__hi.x, v.__lo.x, w.__hi.x, w.__lo.x);
    __dssub(v.__hi.y, v.__lo.y, v.__hi.y, v.__lo.y, w.__hi.y, w.__lo.y);
    __dssub(v.__hi.z, v.__lo.z, v.__hi.z, v.__lo.z, w.__hi.z, w.__lo.z);
    __dssub(v.__hi.w, v.__lo.w, v.__hi.w, v.__lo.w, w.__hi.w, w.__lo.w);
    return v;
}

/**
 * assignment by scalar multiplication
 */
__device__ inline vector<dfloat, 2>& operator*=(vector<dfloat, 2>& v, dfloat const& s)
{
    __dsmul(v.__hi.x, v.__lo.x, v.__hi.x, v.__lo.x, s.__hi, s.__lo);
    __dsmul(v.__hi.y, v.__lo.y, v.__hi.y, v.__lo.y, s.__hi, s.__lo);
    return v;
}

__device__ inline vector<dfloat, 3>& operator*=(vector<dfloat, 3>& v, dfloat const& s)
{
    __dsmul(v.__hi.x, v.__lo.x, v.__hi.x, v.__lo.x, s.__hi, s.__lo);
    __dsmul(v.__hi.y, v.__lo.y, v.__hi.y, v.__lo.y, s.__hi, s.__lo);
    __dsmul(v.__hi.z, v.__lo.z, v.__hi.z, v.__lo.z, s.__hi, s.__lo);
    return v;
}

__device__ inline vector<dfloat, 4>& operator*=(vector<dfloat, 4>& v, dfloat const& s)
{
    __dsmul(v.__hi.x, v.__lo.x, v.__hi.x, v.__lo.x, s.__hi, s.__lo);
    __dsmul(v.__hi.y, v.__lo.y, v.__hi.y, v.__lo.y, s.__hi, s.__lo);
    __dsmul(v.__hi.z, v.__lo.z, v.__hi.z, v.__lo.z, s.__hi, s.__lo);
    __dsmul(v.__hi.w, v.__lo.w, v.__hi.w, v.__lo.w, s.__hi, s.__lo);
    return v;
}

/**
 * assignment by scalar division
 */
__device__ inline vector<dfloat, 2>& operator/=(vector<dfloat, 2>& v, dfloat const& s)
{
    __dsdiv(v.__hi.x, v.__lo.x, v.__hi.x, v.__lo.x, s.__hi, s.__lo);
    __dsdiv(v.__hi.y, v.__lo.y, v.__hi.y, v.__lo.y, s.__hi, s.__lo);
    return v;
}

__device__ inline vector<dfloat, 3>& operator/=(vector<dfloat, 3>& v, dfloat const& s)
{
    __dsdiv(v.__hi.x, v.__lo.x, v.__hi.x, v.__lo.x, s.__hi, s.__lo);
    __dsdiv(v.__hi.y, v.__lo.y, v.__hi.y, v.__lo.y, s.__hi, s.__lo);
    __dsdiv(v.__hi.z, v.__lo.z, v.__hi.z, v.__lo.z, s.__hi, s.__lo);
    return v;
}

__device__ inline vector<dfloat, 4>& operator/=(vector<dfloat, 4>& v, dfloat const& s)
{
    __dsdiv(v.__hi.x, v.__lo.x, v.__hi.x, v.__lo.x, s.__hi, s.__lo);
    __dsdiv(v.__hi.y, v.__lo.y, v.__hi.y, v.__lo.y, s.__hi, s.__lo);
    __dsdiv(v.__hi.z, v.__lo.z, v.__hi.z, v.__lo.z, s.__hi, s.__lo);
    __dsdiv(v.__hi.w, v.__lo.w, v.__hi.w, v.__lo.w, s.__hi, s.__lo);
    return v;
}

/**
 * scalar product
 */
__device__ inline dfloat operator*(vector<dfloat, 2> v, vector<dfloat, 2> const& w)
{
    __dsmul(v.__hi.x, v.__lo.x, v.__hi.x, v.__lo.x, w.__hi.x, w.__lo.x);
    __dsmul(v.__hi.y, v.__lo.y, v.__hi.y, v.__lo.y, w.__hi.y, w.__lo.y);
    __dsadd(v.__hi.x, v.__lo.x, v.__hi.x, v.__lo.x, v.__hi.y, v.__lo.y);
    return dfloat(v.__hi.x, v.__lo.x);
}

__device__ inline dfloat operator*(vector<dfloat, 3> v, vector<dfloat, 3> const& w)
{
    __dsmul(v.__hi.x, v.__lo.x, v.__hi.x, v.__lo.x, w.__hi.x, w.__lo.x);
    __dsmul(v.__hi.y, v.__lo.y, v.__hi.y, v.__lo.y, w.__hi.y, w.__lo.y);
    __dsmul(v.__hi.z, v.__lo.z, v.__hi.z, v.__lo.z, w.__hi.z, w.__lo.z);
    __dsadd(v.__hi.x, v.__lo.x, v.__hi.x, v.__lo.x, v.__hi.y, v.__lo.y);
    __dsadd(v.__hi.x, v.__lo.x, v.__hi.x, v.__lo.x, v.__hi.z, v.__lo.z);
    return dfloat(v.__hi.x, v.__lo.x);
}

__device__ inline dfloat operator*(vector<dfloat, 4> v, vector<dfloat, 4> const& w)
{
    __dsmul(v.__hi.x, v.__lo.x, v.__hi.x, v.__lo.x, w.__hi.x, w.__lo.x);
    __dsmul(v.__hi.y, v.__lo.y, v.__hi.y, v.__lo.y, w.__hi.y, w.__lo.y);
    __dsmul(v.__hi.z, v.__lo.z, v.__hi.z, v.__lo.z, w.__hi.z, w.__lo.z);
    __dsmul(v.__hi.w, v.__lo.w, v.__hi.w, v.__lo.w, w.__hi.w, w.__lo.w);
    __dsadd(v.__hi.x, v.__lo.x, v.__hi.x, v.__lo.x, v.__hi.y, v.__lo.y);
    __dsadd(v.__hi.x, v.__lo.x, v.__hi.x, v.__lo.x, v.__hi.z, v.__lo.z);
    __dsadd(v.__hi.x, v.__lo.x, v.__hi.x, v.__lo.x, v.__hi.w, v.__lo.w);
    return dfloat(v.__hi.x, v.__lo.x);
}

/**
 * componentwise vector addition
 */
template <unsigned int dimension>
__device__ inline vector<dfloat, dimension> operator+(vector<dfloat, dimension> v, vector<dfloat, dimension> const& w)
{
    v += w;
    return v;
}

/**
 * componentwise vector subtraction
 */
template <unsigned int dimension>
__device__ inline vector<dfloat, dimension> operator-(vector<dfloat, dimension> v, vector<dfloat, dimension> const& w)
{
    v -= w;
    return v;
}

/**
 * scalar multiplication
 */
template <unsigned int dimension>
__device__ inline vector<dfloat, dimension> operator*(vector<dfloat, dimension> v, dfloat const& s)
{
    v *= s;
    return v;
}

/**
 * scalar multiplication
 */
template <unsigned int dimension>
__device__ inline vector<dfloat, dimension> operator*(dfloat const& s, vector<dfloat, dimension> v)
{
    v *= s;
    return v;
}

}} // namespace ljgpu::cu

#endif /* ! LJGPU_MATH_GPU_DSVECTOR_CUH */
