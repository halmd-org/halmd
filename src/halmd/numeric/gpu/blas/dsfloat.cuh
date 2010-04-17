/*
 * Copyright © 2008-2010  Peter Colberg
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

#ifndef HALMD_NUMERIC_GPU_BLAS_DSFLOAT_CUH
#define HALMD_NUMERIC_GPU_BLAS_DSFLOAT_CUH

#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>

#include <halmd/numeric/gpu/blas/detail/dsfun.cuh>

namespace halmd { namespace numeric { namespace gpu { namespace blas
{

/**
 * Double-single floating point value
 */
struct dsfloat
{
    float hi, lo;

    __device__ __host__ dsfloat()
    {}

    __device__ __host__ dsfloat(float a0, float a1)
      : hi(a0), lo(a1)
    {}

    template <typename T>
    __device__ __host__ dsfloat(T a0,
      typename boost::enable_if<boost::is_same<T, float> >::type* dummy = 0)
    {
        detail::dsfeq(hi, lo, a0);
    }

    template <typename T>
    __device__ __host__ dsfloat(T a0,
      typename boost::enable_if<boost::is_integral<T> >::type* dummy = 0)
    {
        detail::dsfeq(hi, lo, a0);
    }

    template <typename T>
    __device__ __host__ dsfloat(T a,
      typename boost::enable_if<boost::is_same<T, double> >::type* dummy = 0)
    {
        detail::dsdeq(hi, lo, a);
    }

    template <typename T>
    __device__ __host__ operator T() const
    {
        return hi;
    }

    __device__ __host__ operator double() const
    {
        return static_cast<double>(hi) + lo;
    }
};

#ifdef __CUDACC__

/**
 * Returns "high" single precision floating-point value
 */
__device__ inline float dsfloat_hi(dsfloat const& v)
{
    return v.hi;
}

/**
 * Returns "low" single precision floating-point value
 */
__device__ inline float dsfloat_lo(dsfloat const& v)
{
    return v.lo;
}

/**
 * Addition by assignment
 */
__device__ inline dsfloat& operator+=(dsfloat& v, dsfloat const& w)
{
    detail::dsadd(v.hi, v.lo, v.hi, v.lo, w.hi, w.lo);
    return v;
}

/**
 * Subtraction by assignment
 */
__device__ inline dsfloat& operator-=(dsfloat& v, dsfloat const& w)
{
    detail::dssub(v.hi, v.lo, v.hi, v.lo, w.hi, w.lo);
    return v;
}

/**
 * Multiplication by assignment
 */
__device__ inline dsfloat& operator*=(dsfloat& v, dsfloat const& w)
{
    detail::dsmul(v.hi, v.lo, v.hi, v.lo, w.hi, w.lo);
    return v;
}

/**
 * Division by assignment
 */
__device__ inline dsfloat& operator/=(dsfloat& v, dsfloat const& w)
{
    detail::dsdiv(v.hi, v.lo, v.hi, v.lo, w.hi, w.lo);
    return v;
}

/**
 * Addition
 */
__device__ inline dsfloat operator+(dsfloat v, dsfloat const& w)
{
    v += w;
    return v;
}

template <typename T>
__device__ inline typename boost::enable_if<boost::is_arithmetic<T>, dsfloat>::type
operator+(T const& v, dsfloat const& w)
{
    return static_cast<dsfloat>(v) + w;
}

template <typename T>
__device__ inline typename boost::enable_if<boost::is_arithmetic<T>, dsfloat>::type
operator+(dsfloat const& v, T const& w)
{
    return v + static_cast<dsfloat>(w);
}

/**
 * Subtraction
 */
__device__ inline dsfloat operator-(dsfloat v, dsfloat const& w)
{
    v -= w;
    return v;
}

template <typename T>
__device__ inline typename boost::enable_if<boost::is_arithmetic<T>, dsfloat>::type
operator-(T const& v, dsfloat const& w)
{
    return static_cast<dsfloat>(v) - w;
}

template <typename T>
__device__ inline typename boost::enable_if<boost::is_arithmetic<T>, dsfloat>::type
operator-(dsfloat const& v, T const& w)
{
    return v - static_cast<dsfloat>(w);
}

/**
 * Multiplication
 */
__device__ inline dsfloat operator*(dsfloat v, dsfloat const& w)
{
    v *= w;
    return v;
}

template <typename T>
__device__ inline typename boost::enable_if<boost::is_arithmetic<T>, dsfloat>::type
operator*(T const& v, dsfloat const& w)
{
    return static_cast<dsfloat>(v) * w;
}

template <typename T>
__device__ inline typename boost::enable_if<boost::is_arithmetic<T>, dsfloat>::type
operator*(dsfloat const& v, T const& w)
{
    return v * static_cast<dsfloat>(w);
}

__device__ inline dsfloat operator/(dsfloat v, dsfloat const& w)
{
    v /= w;
    return v;
}

/**
 * Division
 */
template <typename T>
__device__ inline typename boost::enable_if<boost::is_arithmetic<T>, dsfloat>::type
operator/(T const& v, dsfloat const& w)
{
    return static_cast<dsfloat>(v) / w;
}

template <typename T>
__device__ inline typename boost::enable_if<boost::is_arithmetic<T>, dsfloat>::type
operator/(dsfloat const& v, T const& w)
{
    return v / static_cast<dsfloat>(w);
}

/**
 * Square root function
 */
__device__ inline dsfloat sqrt(dsfloat v)
{
    dsfloat w;
    detail::dssqrt(w.hi, w.lo, v.hi, v.lo);
    return w;
}

#endif /* __CUDACC__ */

}}}} // namespace halmd::numeric::gpu::blas

#endif /* ! HALMD_NUMERIC_GPU_BLAS_DSFLOAT_CUH */
