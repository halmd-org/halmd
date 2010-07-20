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

#ifndef HALMD_NUMERIC_GPU_BLAS_DETAIL_VECTOR4D_CUH
#define HALMD_NUMERIC_GPU_BLAS_DETAIL_VECTOR4D_CUH

#include <boost/mpl/and.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/utility/enable_if.hpp>
#include <cuda_runtime.h>

#include <halmd/numeric/gpu/blas/detail/dsfloat.cuh>
#include <halmd/numeric/gpu/blas/detail/storage.cuh>

namespace halmd
{
namespace numeric { namespace gpu { namespace blas
{
namespace detail
{

template <typename T, size_t N>
struct vector;

/**
 * Four-dimensional vector of arbitrary value type
 */
template <typename T>
struct vector<T, 4> : bounded_array<T, 4>
{
    typedef bounded_array<T, 4> _Base;
    typedef typename _Base::value_type value_type;
    enum { static_size = _Base::static_size };

    __device__ vector()
    {}

    /**
     * Initialization by scalar
     */
    template <typename T_>
    __device__ vector(T_ const& s,
      typename boost::enable_if<boost::is_convertible<T_, T> >::type* dummy = 0)
    {
        (*this)[0] = s;
        (*this)[1] = s;
        (*this)[2] = s;
        (*this)[3] = s;
    }

    /**
     * Explicit conversion from vector of convertible element type
     */
    template <typename T_>
    __device__ explicit vector(vector<T_, 4> const& v,
      typename boost::enable_if<boost::is_convertible<T_, T> >::type* dummy = 0)
    {
        (*this)[0] = static_cast<T>(v[0]);
        (*this)[1] = static_cast<T>(v[1]);
        (*this)[2] = static_cast<T>(v[2]);
        (*this)[3] = static_cast<T>(v[3]);
    }
};

/**
 * Four-dimensional single precision floating-point vector
 */
template <>
struct vector<float, 4> : bounded_array<float, 4>
{
    typedef bounded_array<float, 4> _Base;
    typedef _Base::value_type value_type;
    enum { static_size = _Base::static_size };

    __device__ vector() {}

    /**
     * Initialization by scalar
     */
    __device__ vector(float const& s)
    {
        (*this)[0] = s;
        (*this)[1] = s;
        (*this)[2] = s;
        (*this)[3] = s;
    }

    /**
     * Explicit conversion from vector of convertible element type
     */
    template <typename T_>
    __device__ explicit vector(vector<T_, 4> const& v,
      typename boost::enable_if<boost::is_convertible<T_, float> >::type* dummy = 0)
    {
        (*this)[0] = static_cast<float>(v[0]);
        (*this)[1] = static_cast<float>(v[1]);
        (*this)[2] = static_cast<float>(v[2]);
        (*this)[3] = static_cast<float>(v[3]);
    }

    /**
     * Convert from CUDA vector type
     */
    __device__ vector(float4 const& v)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
        (*this)[2] = v.z;
        (*this)[3] = v.w;
    }

    /**
     * Convert to CUDA vector type
     */
    __device__ operator float4() const
    {
        float4 v;
        v.x = (*this)[0];
        v.y = (*this)[1];
        v.z = (*this)[2];
        v.w = (*this)[3];
        return v;
    }
};

/**
 * Four-dimensional unsigned integer vector
 */
template <>
struct vector<unsigned int, 4> : bounded_array<unsigned int, 4>
{
    typedef bounded_array<unsigned int, 4> _Base;
    typedef _Base::value_type value_type;
    enum { static_size = _Base::static_size };

    __device__ vector() {}

    /**
     * Initialization by scalar
     */
    __device__ vector(unsigned int const& s)
    {
        (*this)[0] = s;
        (*this)[1] = s;
        (*this)[2] = s;
        (*this)[3] = s;
    }

    /**
     * Explicit conversion from vector of convertible element type
     */
    template <typename T_>
    __device__ explicit vector(vector<T_, 4> const& v,
      typename boost::enable_if<boost::is_convertible<T_, unsigned int> >::type* dummy = 0)
    {
        (*this)[0] = static_cast<unsigned int>(v[0]);
        (*this)[1] = static_cast<unsigned int>(v[1]);
        (*this)[2] = static_cast<unsigned int>(v[2]);
        (*this)[3] = static_cast<unsigned int>(v[3]);
    }

    /**
     * Convert from CUDA vector type
     */
    __device__ vector(uint4 const& v)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
        (*this)[2] = v.z;
        (*this)[3] = v.w;
    }

    /**
     * Convert to CUDA vector type
     */
    __device__ operator uint4() const
    {
        uint4 v;
        v.x = (*this)[0];
        v.y = (*this)[1];
        v.z = (*this)[2];
        v.w = (*this)[3];
        return v;
    }
};

/**
 * Four-dimensional integer vector
 */
template <>
struct vector<int, 4> : bounded_array<int, 4>
{
    typedef bounded_array<int, 4> _Base;
    typedef _Base::value_type value_type;
    enum { static_size = _Base::static_size };

    __device__ vector() {}

    /**
     * Initialization by scalar
     */
    __device__ vector(int const& s)
    {
        (*this)[0] = s;
        (*this)[1] = s;
        (*this)[2] = s;
        (*this)[3] = s;
    }

    /**
     * Explicit conversion from vector of convertible element type
     */
    template <typename T_>
    __device__ explicit vector(vector<T_, 4> const& v,
      typename boost::enable_if<boost::is_convertible<T_, int> >::type* dummy = 0)
    {
        (*this)[0] = static_cast<int>(v[0]);
        (*this)[1] = static_cast<int>(v[1]);
        (*this)[2] = static_cast<int>(v[2]);
        (*this)[3] = static_cast<int>(v[3]);
    }

    /**
     * Convert from CUDA vector type
     */
    __device__ vector(int4 const& v)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
        (*this)[2] = v.z;
        (*this)[3] = v.w;
    }

    /**
     * Convert to CUDA vector type
     */
    __device__ operator int4() const
    {
        int4 v;
        v.x = (*this)[0];
        v.y = (*this)[1];
        v.z = (*this)[2];
        v.w = (*this)[3];
        return v;
    }
};

/**
 * Four-dimensional double-single precision floating-point vector
 */
template <>
struct vector<dsfloat, 4> : bounded_array<dsfloat, 4>
{
    typedef bounded_array<dsfloat, 4> _Base;
    typedef _Base::value_type value_type;
    enum { static_size = _Base::static_size };

    __device__ vector() {}

    /**
     * Initialization by scalar
     */
    template <typename T_>
    __device__ vector(T_ const& s,
      typename boost::enable_if<boost::is_convertible<T_, dsfloat> >::type* dummy = 0)
    {
        (*this)[0] = s;
        (*this)[1] = s;
        (*this)[2] = s;
        (*this)[3] = s;
    }

    /**
     * Implicit conversion from vector of convertible element type
     */
    template <typename T_>
    __device__ vector(vector<T_, 4> const& v,
      typename boost::enable_if<boost::is_convertible<T_, dsfloat> >::type* dummy = 0)
    {
        (*this)[0] = static_cast<dsfloat>(v[0]);
        (*this)[1] = static_cast<dsfloat>(v[1]);
        (*this)[2] = static_cast<dsfloat>(v[2]);
        (*this)[3] = static_cast<dsfloat>(v[3]);
    }

    __device__ vector(vector<float, 4> const& v, vector<float, 4> const& w)
    {
        (*this)[0] = dsfloat(v[0], w[0]);
        (*this)[1] = dsfloat(v[1], w[1]);
        (*this)[2] = dsfloat(v[2], w[2]);
        (*this)[3] = dsfloat(v[3], w[3]);
    }
};

#ifdef __CUDACC__

/**
 * Returns "high" and "low" single precision vector tuple
 */
__device__ inline tuple<vector<float, 4>, vector<float, 4> > split(vector<dsfloat, 4> const& v)
{
    vector<float, 4> hi, lo;
    tie(hi[0], lo[0]) = split(v[0]);
    tie(hi[1], lo[1]) = split(v[1]);
    tie(hi[2], lo[2]) = split(v[2]);
    tie(hi[3], lo[3]) = split(v[3]);
    return make_tuple(hi, lo);
}

/**
 * Assignment by elementwise vector addition
 */
template <typename T, typename T_>
__device__ inline typename boost::enable_if<boost::is_convertible<T_, T>, vector<T, 4>&>::type
operator+=(vector<T, 4>& v, vector<T_, 4> const& w)
{
    v[0] += w[0];
    v[1] += w[1];
    v[2] += w[2];
    v[3] += w[3];
    return v;
}

/**
 * Assignment by elementwise vector subtraction
 */
template <typename T, typename T_>
__device__ inline typename boost::enable_if<boost::is_convertible<T_, T>, vector<T, 4>&>::type
operator-=(vector<T, 4>& v, vector<T_, 4> const& w)
{
    v[0] -= w[0];
    v[1] -= w[1];
    v[2] -= w[2];
    v[3] -= w[3];
    return v;
}

/**
 * Assignment by scalar multiplication
 */
template <typename T, typename T_>
__device__ inline typename boost::enable_if<boost::is_convertible<T_, T>, vector<T, 4>&>::type
operator*=(vector<T, 4>& v, T_ s)
{
    v[0] *= s;
    v[1] *= s;
    v[2] *= s;
    v[3] *= s;
    return v;
}

/**
 * Assignment by scalar division
 */
template <typename T, typename T_>
__device__ inline typename boost::enable_if<boost::is_convertible<T_, T>, vector<T, 4>&>::type
operator/=(vector<T, 4>& v, T_ s)
{
    v[0] /= s;
    v[1] /= s;
    v[2] /= s;
    v[3] /= s;
    return v;
}

/**
 * Assignment by scalar modulus
 */
template <typename T, typename T_>
__device__ inline typename boost::enable_if<boost::mpl::and_<boost::is_integral<T>, boost::is_integral<T_> >, vector<T, 4>&>::type
operator%=(vector<T, 4>& v, T_ s)
{
    v[0] %= s;
    v[1] %= s;
    v[2] %= s;
    v[3] %= s;
    return v;
}

/**
 * Elementwise vector addition
 */
template <typename T>
__device__ inline vector<T, 4> operator+(vector<T, 4> v, vector<T, 4> const& w)
{
    v += w;
    return v;
}

/**
 * Elementwise vector subtraction
 */
template <typename T>
__device__ inline vector<T, 4> operator-(vector<T, 4> v, vector<T, 4> const& w)
{
    v -= w;
    return v;
}

/**
 * Elementwise change of sign
 */
template <typename T>
__device__ inline vector<T, 4> operator-(vector<T, 4> v)
{
    v[0] = -v[0];
    v[1] = -v[1];
    v[2] = -v[2];
    v[3] = -v[3];
    return v;
}

/**
 * Scalar multiplication
 */
template <typename T, typename T_>
__device__ inline typename boost::enable_if<boost::is_convertible<T_, T>, vector<T, 4> >::type
operator*(vector<T, 4> v, T_ s)
{
    v[0] *= s;
    v[1] *= s;
    v[2] *= s;
    v[3] *= s;
    return v;
}

/**
 * Scalar multiplication
 */
template <typename T, typename T_>
__device__ inline typename boost::enable_if<boost::is_convertible<T_, T>, vector<T, 4> >::type
operator*(T_ s, vector<T, 4> v)
{
    v[0] *= s;
    v[1] *= s;
    v[2] *= s;
    v[3] *= s;
    return v;
}

/**
 * Scalar division
 */
template <typename T, typename T_>
__device__ inline typename boost::enable_if<boost::is_convertible<T_, T>, vector<T, 4> >::type
operator/(vector<T, 4> v, T_ s)
{
    v[0] /= s;
    v[1] /= s;
    v[2] /= s;
    v[3] /= s;
    return v;
}

/**
 * Scalar modulus
 */
template <typename T, typename T_>
__device__ inline typename boost::enable_if<boost::mpl::and_<boost::is_integral<T>, boost::is_integral<T_> >, vector<T, 4> >::type
operator%(vector<T, 4> v, T_ s)
{
    v[0] %= s;
    v[1] %= s;
    v[2] %= s;
    v[3] %= s;
    return v;
}

/**
 * Inner product
 */
template <typename T>
__device__ inline T inner_prod(vector<T, 4> const& v, vector<T, 4> const& w)
{
    T s = v[0] * w[0];
    s  += v[1] * w[1];
    s  += v[2] * w[2];
    s  += v[3] * w[3];
    return s;
}

/**
 * Elementwise vector multiplication
 */
template <typename T>
__device__ inline vector<T, 4> element_prod(vector<T, 4> v, vector<T, 4> const& w)
{
    v[0] *= w[0];
    v[1] *= w[1];
    v[2] *= w[2];
    v[3] *= w[3];
    return v;
}

/**
 * Elementwise vector division
 */
template <typename T>
__device__ inline vector<T, 4> element_div(vector<T, 4> v, vector<T, 4> const& w)
{
    v[0] /= w[0];
    v[1] /= w[1];
    v[2] /= w[2];
    v[3] /= w[3];
    return v;
}

/**
 * Elementwise vector modulus
 */
template <typename T>
__device__ inline typename boost::enable_if<boost::is_integral<T>, vector<T, 4> >::type
element_mod(vector<T, 4> v, vector<T, 4> const& w)
{
    v[0] %= w[0];
    v[1] %= w[1];
    v[2] %= w[2];
    v[3] %= w[3];
    return v;
}

/**
 * Elementwise round to nearest integer
 */
__device__ inline vector<float, 4> rint(vector<float, 4> v)
{
    v[0] = ::rintf(v[0]);
    v[1] = ::rintf(v[1]);
    v[2] = ::rintf(v[2]);
    v[3] = ::rintf(v[3]);
    return v;
}

/**
 * Elementwise round to nearest integer, away from zero
 */
__device__ inline vector<float, 4> round(vector<float, 4> v)
{
    v[0] = ::roundf(v[0]);
    v[1] = ::roundf(v[1]);
    v[2] = ::roundf(v[2]);
    v[3] = ::roundf(v[3]);
    return v;
}

/**
 * Elementwise round to nearest integer not greater than argument
 */
__device__ inline vector<float, 4> floor(vector<float, 4> v)
{
    v[0] = ::floorf(v[0]);
    v[1] = ::floorf(v[1]);
    v[2] = ::floorf(v[2]);
    v[3] = ::floorf(v[3]);
    return v;
}

/**
 * Elementwise round to nearest integer not less argument
 */
__device__ inline vector<float, 4> ceil(vector<float, 4> v)
{
    v[0] = ::ceilf(v[0]);
    v[1] = ::ceilf(v[1]);
    v[2] = ::ceilf(v[2]);
    v[3] = ::ceilf(v[3]);
    return v;
}

/**
 * Elementwise square root function
 */
__device__ inline vector<float, 4> sqrt(vector<float, 4> v)
{
    v[0] = ::sqrtf(v[0]);
    v[1] = ::sqrtf(v[1]);
    v[2] = ::sqrtf(v[2]);
    v[3] = ::sqrtf(v[3]);
    return v;
}

/**
 * Elementwise cosine function
 */
__device__ inline vector<float, 4> cos(vector<float, 4> v)
{
    v[0] = ::cosf(v[0]);
    v[1] = ::cosf(v[1]);
    v[2] = ::cosf(v[2]);
    v[3] = ::cosf(v[3]);
    return v;
}

/**
 * Elementwise sine function
 */
__device__ inline vector<float, 4> sin(vector<float, 4> v)
{
    v[0] = ::sinf(v[0]);
    v[1] = ::sinf(v[1]);
    v[2] = ::sinf(v[2]);
    v[3] = ::sinf(v[3]);
    return v;
}

/**
 * Elementwise absolute value
 */
__device__ inline vector<float, 4> fabs(vector<float, 4> v)
{
    v[0] = ::fabsf(v[0]);
    v[1] = ::fabsf(v[1]);
    v[2] = ::fabsf(v[2]);
    v[3] = ::fabsf(v[3]);
    return v;
}

/**
 * Convert floating-point elements to integers, rounding to nearest even integer
 */
__device__ inline vector<int, 4> __float2int_rn(vector<float, 4> const& v)
{
    vector<int, 4> w;
    w[0] = ::__float2int_rn(v[0]);
    w[1] = ::__float2int_rn(v[1]);
    w[2] = ::__float2int_rn(v[2]);
    w[3] = ::__float2int_rn(v[3]);
    return w;
}

/**
 * Convert floating-point elements to integers, rounding towards zero
 */
__device__ inline vector<int, 4> __float2int_rz(vector<float, 4> const& v)
{
    vector<int, 4> w;
    w[0] = ::__float2int_rz(v[0]);
    w[1] = ::__float2int_rz(v[1]);
    w[2] = ::__float2int_rz(v[2]);
    w[3] = ::__float2int_rz(v[3]);
    return w;
}

/**
 * Convert floating-point elements to integers, rounding to positive infinity
 */
__device__ inline vector<int, 4> __float2int_ru(vector<float, 4> const& v)
{
    vector<int, 4> w;
    w[0] = ::__float2int_ru(v[0]);
    w[1] = ::__float2int_ru(v[1]);
    w[2] = ::__float2int_ru(v[2]);
    w[3] = ::__float2int_ru(v[3]);
    return w;
}

/**
 * Convert floating-point elements to integers, rounding to negative infinity
 */
__device__ inline vector<int, 4> __float2int_rd(vector<float, 4> const& v)
{
    vector<int, 4> w;
    w[0] = ::__float2int_rd(v[0]);
    w[1] = ::__float2int_rd(v[1]);
    w[2] = ::__float2int_rd(v[2]);
    w[3] = ::__float2int_rd(v[3]);
    return w;
}

/**
 * Convert floating-point elements to unsigned integers, rounding to nearest even integer
 */
__device__ inline vector<unsigned int, 4> __float2uint_rn(vector<float, 4> const& v)
{
    vector<unsigned int, 4> w;
    w[0] = ::__float2uint_rn(v[0]);
    w[1] = ::__float2uint_rn(v[1]);
    w[2] = ::__float2uint_rn(v[2]);
    w[3] = ::__float2uint_rn(v[3]);
    return w;
}

/**
 * Convert floating-point elements to unsigned integers, rounding towards zero
 */
__device__ inline vector<unsigned int, 4> __float2uint_rz(vector<float, 4> const& v)
{
    vector<unsigned int, 4> w;
    w[0] = ::__float2uint_rz(v[0]);
    w[1] = ::__float2uint_rz(v[1]);
    w[2] = ::__float2uint_rz(v[2]);
    w[3] = ::__float2uint_rz(v[3]);
    return w;
}

/**
 * Convert floating-point elements to unsigned integers, rounding to positive infinity
 */
__device__ inline vector<unsigned int, 4> __float2uint_ru(vector<float, 4> const& v)
{
    vector<unsigned int, 4> w;
    w[0] = ::__float2uint_ru(v[0]);
    w[1] = ::__float2uint_ru(v[1]);
    w[2] = ::__float2uint_ru(v[2]);
    w[3] = ::__float2uint_ru(v[3]);
    return w;
}

/**
 * Convert floating-point elements to unsigned integers, rounding to negative infinity
 */
__device__ inline vector<unsigned int, 4> __float2uint_rd(vector<float, 4> const& v)
{
    vector<unsigned int, 4> w;
    w[0] = ::__float2uint_rd(v[0]);
    w[1] = ::__float2uint_rd(v[1]);
    w[2] = ::__float2uint_rd(v[2]);
    w[3] = ::__float2uint_rd(v[3]);
    return w;
}

/**
 * Limit floating-point elements to unit interval [0, 1]
 */
__device__ inline vector<float, 4> __saturate(vector<float, 4> v)
{
    v[0] = ::__saturatef(v[0]);
    v[1] = ::__saturatef(v[1]);
    v[2] = ::__saturatef(v[2]);
    v[3] = ::__saturatef(v[3]);
    return v;
}

/**
 * Floating-point remainder function, round towards nearest integer
 */
__device__ inline vector<float, 4> remainder(vector<float, 4> v, vector<float, 4> const& w)
{
    v[0] = ::remainderf(v[0], w[0]);
    v[1] = ::remainderf(v[1], w[1]);
    v[2] = ::remainderf(v[2], w[2]);
    v[3] = ::remainderf(v[3], w[3]);
    return v;
}

/**
 * Floating-point remainder function, round towards zero
 */
__device__ inline vector<float, 4> fmod(vector<float, 4> v, vector<float, 4> const& w)
{
    v[0] = ::fmodf(v[0], w[0]);
    v[1] = ::fmodf(v[1], w[1]);
    v[2] = ::fmodf(v[2], w[2]);
    v[3] = ::fmodf(v[3], w[3]);
    return v;
}

/**
 * Fast, accurate floating-point division by s < 2^126
 */
__device__ inline vector<float, 4> __fdivide(vector<float, 4> v, vector<float, 4> const& w)
{
    v[0] = ::__fdividef(v[0], w[0]);
    v[1] = ::__fdividef(v[1], w[1]);
    v[2] = ::__fdividef(v[2], w[2]);
    v[3] = ::__fdividef(v[3], w[3]);
    return v;
}

#endif /* __CUDACC__ */

} // namespace detail

}}} // namespace numeric::gpu::blas

} // namespace halmd

#endif /* ! HALMD_NUMERIC_GPU_BLAS_DETAIL_VECTOR4D_CUH */
