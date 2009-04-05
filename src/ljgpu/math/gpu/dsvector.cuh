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

template <unsigned int dimension>
struct vector<dfloat, dimension>
{
    enum { static_size = dimension };
    typedef vector<float, dimension> _Base;
    typedef dfloat value_type;

    _Base f0, f1;

    __device__ inline vector() {}

    __device__ inline vector(dfloat s) : f0(s.f0), f1(s.f1) {}

    __device__ inline vector(float s) : f0(s), f1(0) {}

    __device__ inline vector(_Base v, _Base w) : f0(v), f1(w) {}

    __device__ inline vector(_Base v) : f0(v), f1(0) {}

    __device__ inline operator _Base() const { return f0; }
};

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
