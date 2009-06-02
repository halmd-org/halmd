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
    float2 f0, f1;

    __device__ inline vector() {}

    __device__ inline vector(dfloat s) : f0(make_float2(s.f0, s.f0)), f1(make_float2(s.f1, s.f1)) {}

    __device__ inline vector(float s) : f0(make_float2(s, s)), f1(make_float2(0, 0)) {}

    __device__ inline vector(vector<float, 2> v, vector<float, 2> w) : f0(v), f1(w) {}

    __device__ inline vector(vector<float, 2> v) : f0(v), f1(make_float2(0, 0)) {}

    __device__ inline operator vector<float, 2>() const { return f0; }
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
    float3 f0, f1;

    __device__ inline vector() {}

    __device__ inline vector(dfloat s) : f0(make_float3(s.f0, s.f0, s.f0)), f1(make_float3(s.f1, s.f1, s.f1)) {}

    __device__ inline vector(float s) : f0(make_float3(s, s, s)), f1(make_float3(0, 0, 0)) {}

    __device__ inline vector(vector<float, 3> v, vector<float, 3> w) : f0(v), f1(w) {}

    __device__ inline vector(vector<float, 3> v) : f0(v), f1(make_float3(0, 0, 0)) {}

    __device__ inline operator vector<float, 3>() const { return f0; }
};

/**
 * Double-single precision vector for __shared__ arrays.
 *
 * This class implements a work-around for an NVCC compiler bug:
 * A __shared__ array may only be declared as a struct type if the
 * struct's members are POD or CUDA types, e.g. we cannot use a
 * custom two- or three-dimensional vector of dfloats.
 */
template <typename T, unsigned int dimension>
struct __vector;

template <>
struct __vector<dfloat, 2>
{
    enum { static_size = 2 };

    float2 f0, f1;

    /**
     * Convert from algebraic to non-algebraic vector.
     *
     * This class may not have a non-default constructor because
     * of an NVCC compiler bug, so we define the assignment operator
     * instead.
     */
    __device__ inline __vector<dfloat, 2>& operator=(vector<dfloat, 2> v)
    {
	f0.x = v.f0.x;
	f0.y = v.f0.y;
	f1.x = v.f1.x;
	f1.y = v.f1.y;
	return *this;
    }

    /**
     * Convert from non-algebraic to algebraic vector.
     */
    __device__ inline operator vector<dfloat, 2>()
    {
	vector<dfloat, 2> v;
	v.f0.x = f0.x;
	v.f0.y = f0.y;
	v.f1.x = f1.x;
	v.f1.y = f1.y;
	return v;
    }
};

template <>
struct __vector<dfloat, 3>
{
    enum { static_size = 3 };

    float3 f0, f1;

    /**
     * Convert from algebraic to non-algebraic vector.
     *
     * This class may not have a non-default constructor because
     * of an NVCC compiler bug, so we define the assignment operator
     * instead.
     */
    __device__ inline __vector<dfloat, 3>& operator=(vector<dfloat, 3> v)
    {
	f0.x = v.f0.x;
	f0.y = v.f0.y;
	f0.z = v.f0.z;
	f1.x = v.f1.x;
	f1.y = v.f1.y;
	f1.z = v.f1.z;
	return *this;
    }

    /**
     * Convert from non-algebraic to algebraic vector.
     */
    __device__ inline operator vector<dfloat, 3>()
    {
	vector<dfloat, 3> v;
	v.f0.x = f0.x;
	v.f0.y = f0.y;
	v.f0.z = f0.z;
	v.f1.x = f1.x;
	v.f1.y = f1.y;
	v.f1.z = f1.z;
	return v;
    }
};

/**
 * returns high-word floating point vector
 */
template <unsigned int dimension>
__device__ inline vector<float, dimension> dfloat2hi(vector<dfloat, dimension> const& v)
{
    return v.f0;
}

/**
 * returns low-word floating point vector
 */
template <unsigned int dimension>
__device__ inline vector<float, dimension> dfloat2lo(vector<dfloat, dimension> const& v)
{
    return v.f1;
}

/**
 * assignment by componentwise vector addition
 */
__device__ inline vector<dfloat, 2>& operator+=(vector<dfloat, 2>& v, vector<dfloat, 2> const& w)
{
    __dsadd(v.f0.x, v.f1.x, v.f0.x, v.f1.x, w.f0.x, w.f1.x);
    __dsadd(v.f0.y, v.f1.y, v.f0.y, v.f1.y, w.f0.y, w.f1.y);
    return v;
}

__device__ inline vector<dfloat, 3>& operator+=(vector<dfloat, 3>& v, vector<dfloat, 3> const& w)
{
    __dsadd(v.f0.x, v.f1.x, v.f0.x, v.f1.x, w.f0.x, w.f1.x);
    __dsadd(v.f0.y, v.f1.y, v.f0.y, v.f1.y, w.f0.y, w.f1.y);
    __dsadd(v.f0.z, v.f1.z, v.f0.z, v.f1.z, w.f0.z, w.f1.z);
    return v;
}

/**
 * assignment by componentwise vector subtraction
 */
__device__ inline vector<dfloat, 2>& operator-=(vector<dfloat, 2>& v, vector<dfloat, 2> const& w)
{
    __dssub(v.f0.x, v.f1.x, v.f0.x, v.f1.x, w.f0.x, w.f1.x);
    __dssub(v.f0.y, v.f1.y, v.f0.y, v.f1.y, w.f0.y, w.f1.y);
    return v;
}

__device__ inline vector<dfloat, 3>& operator-=(vector<dfloat, 3>& v, vector<dfloat, 3> const& w)
{
    __dssub(v.f0.x, v.f1.x, v.f0.x, v.f1.x, w.f0.x, w.f1.x);
    __dssub(v.f0.y, v.f1.y, v.f0.y, v.f1.y, w.f0.y, w.f1.y);
    __dssub(v.f0.z, v.f1.z, v.f0.z, v.f1.z, w.f0.z, w.f1.z);
    return v;
}
/**
 * assignment by scalar multiplication
 */
__device__ inline vector<dfloat, 2>& operator*=(vector<dfloat, 2>& v, dfloat const& s)
{
    __dsmul(v.f0.x, v.f1.x, v.f0.x, v.f1.x, s.f0, s.f1);
    __dsmul(v.f0.y, v.f1.y, v.f0.y, v.f1.y, s.f0, s.f1);
    return v;
}

__device__ inline vector<dfloat, 3>& operator*=(vector<dfloat, 3>& v, dfloat const& s)
{
    __dsmul(v.f0.x, v.f1.x, v.f0.x, v.f1.x, s.f0, s.f1);
    __dsmul(v.f0.y, v.f1.y, v.f0.y, v.f1.y, s.f0, s.f1);
    __dsmul(v.f0.z, v.f1.z, v.f0.z, v.f1.z, s.f0, s.f1);
    return v;
}

/**
 * assignment by scalar division
 */
__device__ inline vector<dfloat, 2>& operator/=(vector<dfloat, 2>& v, dfloat const& s)
{
    __dsdiv(v.f0.x, v.f1.x, v.f0.x, v.f1.x, s.f0, s.f1);
    __dsdiv(v.f0.y, v.f1.y, v.f0.y, v.f1.y, s.f0, s.f1);
    return v;
}

__device__ inline vector<dfloat, 3>& operator/=(vector<dfloat, 3>& v, dfloat const& s)
{
    __dsdiv(v.f0.x, v.f1.x, v.f0.x, v.f1.x, s.f0, s.f1);
    __dsdiv(v.f0.y, v.f1.y, v.f0.y, v.f1.y, s.f0, s.f1);
    __dsdiv(v.f0.z, v.f1.z, v.f0.z, v.f1.z, s.f0, s.f1);
    return v;
}

/**
 * scalar product
 */
__device__ inline dfloat operator*(vector<dfloat, 2> v, vector<dfloat, 2> const& w)
{
    __dsmul(v.f0.x, v.f1.x, v.f0.x, v.f1.x, w.f0.x, w.f1.x);
    __dsmul(v.f0.y, v.f1.y, v.f0.y, v.f1.y, w.f0.y, w.f1.y);
    __dsadd(v.f0.x, v.f1.x, v.f0.x, v.f1.x, v.f0.y, v.f1.y);
    return dfloat(v.f0.x, v.f1.x);
}

__device__ inline dfloat operator*(vector<dfloat, 3> v, vector<dfloat, 3> const& w)
{
    __dsmul(v.f0.x, v.f1.x, v.f0.x, v.f1.x, w.f0.x, w.f1.x);
    __dsmul(v.f0.y, v.f1.y, v.f0.y, v.f1.y, w.f0.y, w.f1.y);
    __dsmul(v.f0.z, v.f1.z, v.f0.z, v.f1.z, w.f0.z, w.f1.z);
    __dsadd(v.f0.x, v.f1.x, v.f0.x, v.f1.x, v.f0.y, v.f1.y);
    __dsadd(v.f0.x, v.f1.x, v.f0.x, v.f1.x, v.f0.z, v.f1.z);
    return dfloat(v.f0.x, v.f1.x);
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
