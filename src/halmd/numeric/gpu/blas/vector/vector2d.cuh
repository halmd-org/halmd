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

#ifndef HALMD_NUMERIC_GPU_BLAS_VECTOR_VECTOR2D_CUH
#define HALMD_NUMERIC_GPU_BLAS_VECTOR_VECTOR2D_CUH

#include <boost/type_traits/is_convertible.hpp>
#include <boost/utility/enable_if.hpp>
#include <cuda_runtime.h>

#include <halmd/numeric/gpu/blas/dsfloat.cuh>
#include <halmd/numeric/gpu/blas/vector/storage.cuh>

namespace halmd { namespace numeric { namespace gpu { namespace blas
{

template <typename T, size_t N>
struct vector;

/**
 * Two-dimensional vector of arbitrary value type
 */
template <typename T>
struct vector<T, 2> : bounded_array<T, 2>
{
    typedef bounded_array<T, 2> _Base;
    typedef typename _Base::value_type value_type;
    enum { static_size = _Base::static_size };

    __device__ vector()
    {}

    __device__ vector(T const& s)
    {
        (*this)[0] = s;
        (*this)[1] = s;
    }

    template <typename T_>
    __device__ vector(vector<T_, 2> const& v,
      typename boost::enable_if<boost::is_convertible<T_, T> >::type* dummy = 0)
    {
        (*this)[0] = v[0];
        (*this)[1] = v[1];
    }
};

/**
 * Two-dimensional single precision floating-point vector
 */
template <>
struct vector<float, 2> : bounded_array<float, 2>
{
    typedef bounded_array<float, 2> _Base;
    typedef typename _Base::value_type value_type;
    enum { static_size = _Base::static_size };

    __device__ vector() {}

    __device__ vector(float const& s)
    {
        (*this)[0] = s;
        (*this)[1] = s;
    }

    template <typename T_>
    __device__ vector(vector<T_, 2> const& v,
      typename boost::enable_if<boost::is_convertible<T_, float> >::type* dummy = 0)
    {
        (*this)[0] = v[0];
        (*this)[1] = v[1];
    }

    /**
     * Convert from CUDA vector type
     */
    __device__ vector(float2 const& v)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
    }

    /**
     * Convert to CUDA vector type
     */
    __device__ operator float2() const
    {
        float2 v;
        v.x = (*this)[0];
        v.y = (*this)[1];
        return v;
    }
};

/**
 * Two-dimensional double-single precision floating-point vector
 */
template <>
struct vector<dsfloat, 2> : bounded_array<dsfloat, 2>
{
    typedef bounded_array<dsfloat, 2> _Base;
    typedef typename _Base::value_type value_type;
    enum { static_size = _Base::static_size };

    __device__ vector() {}

    __device__ vector(dsfloat const& s)
    {
        (*this)[0] = s;
        (*this)[1] = s;
    }

    template <typename T_>
    __device__ vector(vector<T_, 2> const& v,
      typename boost::enable_if<boost::is_convertible<T_, dsfloat> >::type* dummy = 0)
    {
        (*this)[0] = v[0];
        (*this)[1] = v[1];
    }

    __device__ vector(vector<float, 2> const& v, vector<float, 2> const& w)
    {
        (*this)[0] = dsfloat(v[0], w[0]);
        (*this)[1] = dsfloat(v[1], w[1]);
    }
};

#ifdef __CUDACC__

/**
 * Returns "high" single precision floating-point vector
 */
__device__ inline vector<float, 2> dsfloat_hi(vector<dsfloat, 2> v)
{
    vector<float, 2> w;
    w[0] = dsfloat_hi(v[0]);
    w[1] = dsfloat_hi(v[1]);
    return w;
}

/**
 * Returns "low" single precision floating-point vector
 */
__device__ inline vector<float, 2> dsfloat_lo(vector<dsfloat, 2> v)
{
    vector<float, 2> w;
    w[0] = dsfloat_lo(v[0]);
    w[1] = dsfloat_lo(v[1]);
    return w;
}

/**
 * Assignment by componentwise vector addition
 */
template <typename T>
// work around NVCC compiler error with __device__ template function that returns reference
__device__ inline typename boost::enable_if<boost::true_type, vector<T, 2>&>::type
operator+=(vector<T, 2>& v, vector<T, 2> const& w)
{
    v[0] += w[0];
    v[1] += w[1];
    return v;
}

/**
 * Assignment by componentwise vector subtraction
 */
template <typename T>
// work around NVCC compiler error with __device__ template function that returns reference
__device__ inline typename boost::enable_if<boost::true_type, vector<T, 2>&>::type
operator-=(vector<T, 2>& v, vector<T, 2> const& w)
{
    v[0] -= w[0];
    v[1] -= w[1];
    return v;
}

/**
 * Assignment by scalar multiplication
 */
template <typename T, typename T_>
__device__ inline typename boost::enable_if<boost::is_convertible<T_, T>, vector<T, 2>&>::type
operator*=(vector<T, 2>& v, T_ s)
{
    v[0] *= s;
    v[1] *= s;
    return v;
}

/**
 * Assignment by scalar division
 */
template <typename T, typename T_>
__device__ inline typename boost::enable_if<boost::is_convertible<T_, T>, vector<T, 2>&>::type
operator/=(vector<T, 2>& v, T_ s)
{
    v[0] /= s;
    v[1] /= s;
    return v;
}

/**
 * Componentwise vector addition
 */
template <typename T>
__device__ inline vector<T, 2> operator+(vector<T, 2> v, vector<T, 2> const& w)
{
    v += w;
    return v;
}

/**
 * Componentwise vector subtraction
 */
template <typename T>
__device__ inline vector<T, 2> operator-(vector<T, 2> v, vector<T, 2> const& w)
{
    v -= w;
    return v;
}

/**
 * Componentwise change of sign
 */
template <typename T>
__device__ inline vector<T, 2> operator-(vector<T, 2> v)
{
    v[0] = -v[0];
    v[1] = -v[1];
    return v;
}

/**
 * Scalar product
 */
template <typename T>
__device__ inline T operator*(vector<T, 2> const& v, vector<T, 2> const& w)
{
    T s = v[0] * w[0];
    s  += v[1] * w[1];
    return s;
}

/**
 * Scalar multiplication
 */
template <typename T, typename T_>
__device__ inline typename boost::enable_if<boost::is_convertible<T_, T>, vector<T, 2> >::type
operator*(vector<T, 2> v, T_ s)
{
    v[0] *= s;
    v[1] *= s;
    return v;
}

/**
 * Scalar multiplication
 */
template <typename T, typename T_>
__device__ inline typename boost::enable_if<boost::is_convertible<T_, T>, vector<T, 2> >::type
operator*(T_ s, vector<T, 2> v)
{
    v[0] *= s;
    v[1] *= s;
    return v;
}

/**
 * Scalar division
 */
template <typename T, typename T_>
__device__ inline typename boost::enable_if<boost::is_convertible<T_, T>, vector<T, 2> >::type
operator/(vector<T, 2> v, T_ s)
{
    v[0] /= s;
    v[1] /= s;
    return v;
}

/**
 * Componentwise round to nearest integer
 */
__device__ inline vector<float, 2> rint(vector<float, 2> v)
{
    v[0] = ::rintf(v[0]);
    v[1] = ::rintf(v[1]);
    return v;
}

/**
 * Componentwise round to nearest integer, away from zero
 */
__device__ vector<float, 2> round(vector<float, 2> v)
{
    v[0] = ::roundf(v[0]);
    v[1] = ::roundf(v[1]);
    return v;
}

/**
 * Componentwise round to nearest integer not greater than argument
 */
__device__ inline vector<float, 2> floor(vector<float, 2> v)
{
    v[0] = ::floorf(v[0]);
    v[1] = ::floorf(v[1]);
    return v;
}

/**
 * Componentwise round to nearest integer not less argument
 */
__device__ inline vector<float, 2> ceil(vector<float, 2> v)
{
    v[0] = ::ceilf(v[0]);
    v[1] = ::ceilf(v[1]);
    return v;
}

/**
 * Componentwise square root function
 */
__device__ inline vector<float, 2> sqrt(vector<float, 2> v)
{
    v[0] = ::sqrtf(v[0]);
    v[1] = ::sqrtf(v[1]);
    return v;
}

/**
 * Componentwise cosine function
 */
__device__ inline vector<float, 2> cos(vector<float, 2> v)
{
    v[0] = ::cosf(v[0]);
    v[1] = ::cosf(v[1]);
    return v;
}

/**
 * Componentwise sine function
 */
__device__ inline vector<float, 2> sin(vector<float, 2> v)
{
    v[0] = ::sinf(v[0]);
    v[1] = ::sinf(v[1]);
    return v;
}

/**
 * Componentwise absolute value
 */
__device__ inline vector<float, 2> fabs(vector<float, 2> v)
{
    v[0] = ::fabsf(v[0]);
    v[1] = ::fabsf(v[1]);
    return v;
}

/**
 * Convert floating-point components to integers, rounding to nearest even integer
 */
__device__ inline vector<int, 2> __float2int_rn(vector<float, 2> const& v)
{
    vector<int, 2> w;
    w[0] = ::__float2int_rn(v[0]);
    w[1] = ::__float2int_rn(v[1]);
    return w;
}

/**
 * Convert floating-point components to integers, rounding towards zero
 */
__device__ inline vector<int, 2> __float2int_rz(vector<float, 2> const& v)
{
    vector<int, 2> w;
    w[0] = ::__float2int_rz(v[0]);
    w[1] = ::__float2int_rz(v[1]);
    return w;
}

/**
 * Convert floating-point components to integers, rounding to positive infinity
 */
__device__ inline vector<int, 2> __float2int_ru(vector<float, 2> const& v)
{
    vector<int, 2> w;
    w[0] = ::__float2int_ru(v[0]);
    w[1] = ::__float2int_ru(v[1]);
    return w;
}

/**
 * Convert floating-point components to integers, rounding to negative infinity
 */
__device__ inline vector<int, 2> __float2int_rd(vector<float, 2> const& v)
{
    vector<int, 2> w;
    w[0] = ::__float2int_rd(v[0]);
    w[1] = ::__float2int_rd(v[1]);
    return w;
}

/**
 * Convert floating-point components to unsigned integers, rounding to nearest even integer
 */
__device__ inline vector<unsigned int, 2> __float2uint_rn(vector<float, 2> const& v)
{
    vector<unsigned int, 2> w;
    w[0] = ::__float2uint_rn(v[0]);
    w[1] = ::__float2uint_rn(v[1]);
    return w;
}

/**
 * Convert floating-point components to unsigned integers, rounding towards zero
 */
__device__ inline vector<unsigned int, 2> __float2uint_rz(vector<float, 2> const& v)
{
    vector<unsigned int, 2> w;
    w[0] = ::__float2uint_rz(v[0]);
    w[1] = ::__float2uint_rz(v[1]);
    return w;
}

/**
 * Convert floating-point components to unsigned integers, rounding to positive infinity
 */
__device__ inline vector<unsigned int, 2> __float2uint_ru(vector<float, 2> const& v)
{
    vector<unsigned int, 2> w;
    w[0] = ::__float2uint_ru(v[0]);
    w[1] = ::__float2uint_ru(v[1]);
    return w;
}

/**
 * Convert floating-point components to unsigned integers, rounding to negative infinity
 */
__device__ inline vector<unsigned int, 2> __float2uint_rd(vector<float, 2> const& v)
{
    vector<unsigned int, 2> w;
    w[0] = ::__float2uint_rd(v[0]);
    w[1] = ::__float2uint_rd(v[1]);
    return w;
}

/**
 * Limit floating-point components to unit interval [0, 1]
 */
__device__ inline vector<float, 2> __saturate(vector<float, 2> v)
{
    v[0] = ::__saturatef(v[0]);
    v[1] = ::__saturatef(v[1]);
    return v;
}

/**
 * Floating-point remainder function, round towards nearest integer
 */
__device__ inline vector<float, 2> remainder(vector<float, 2> v, float s)
{
    v[0] = ::remainderf(v[0], s);
    v[1] = ::remainderf(v[1], s);
    return v;
}

/**
 * Floating-point remainder function, round towards zero
 */
__device__ inline vector<float, 2> fmod(vector<float, 2> v, float s)
{
    v[0] = ::fmodf(v[0], s);
    v[1] = ::fmodf(v[1], s);
    return v;
}

/**
 * Fast, accurate floating-point division by s < 2^126
 */
__device__ inline vector<float, 2> __fdivide(vector<float, 2> v, float s)
{
    v[0] = ::__fdividef(v[0], s);
    v[1] = ::__fdividef(v[1], s);
    return v;
}

#endif /* __CUDACC__ */

}}}} // namespace halmd::numeric::gpu::blas

#endif /* ! HALMD_NUMERIC_GPU_BLAS_VECTOR_VECTOR2D_CUH */
