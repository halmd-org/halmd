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

#ifndef HALMD_NUMERIC_HOST_BLAS_DETAIL_VECTOR4D_HPP
#define HALMD_NUMERIC_HOST_BLAS_DETAIL_VECTOR4D_HPP

#include <boost/array.hpp>
#include <boost/mpl/and.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/utility/enable_if.hpp>
#include <cmath>
#include <iostream>

#ifdef WITH_CUDA
# include <halmd/numeric/gpu/blas/vector.cuh>
#endif

namespace halmd { namespace numeric { namespace host { namespace blas
{

namespace detail
{

template <typename T, size_t N>
struct vector;

/**
 * Four-dimensional vector
 */
template <typename T>
struct vector<T, 4>
  : public boost::array<T, 4>
{
public:
    typedef boost::array<T, 4> _Base;
    typedef typename _Base::value_type value_type;
    enum { static_size = _Base::static_size };

    vector()
    {}

    /**
     * Explicit conversion from vector of convertible element type
     */
    template <typename T_>
    explicit vector(vector<T_, 4> const& v,
      typename boost::enable_if<boost::is_convertible<T_, T> >::type* dummy = 0)
    {
        (*this)[0] = static_cast<T>(v[0]);
        (*this)[1] = static_cast<T>(v[1]);
        (*this)[2] = static_cast<T>(v[2]);
        (*this)[3] = static_cast<T>(v[3]);
    }

    /**
     * Initialization by scalar
     */
    template <typename T_>
    vector(T_ const& s,
      typename boost::enable_if<boost::is_convertible<T_, T> >::type* dummy = 0)
    {
        (*this)[0] = s;
        (*this)[1] = s;
        (*this)[2] = s;
        (*this)[3] = s;
    }

    /**
     * Initialization by scalar elements
     */
    template <typename T_>
    vector(T_ const& x, T_ const& y, T_ const& z,
      typename boost::enable_if<boost::is_convertible<T_, T> >::type* dummy = 0)
    {
        (*this)[0] = x;
        (*this)[1] = y;
        (*this)[2] = z;
        (*this)[3] = z;
    }

#ifdef WITH_CUDA

    /**
     * Convert from GPU vector
     */
    template <typename T_>
    vector(gpu::blas::vector<T_, 4> const& v,
      typename boost::enable_if<boost::is_convertible<T_, T> >::type* dummy = 0)
    {
        (*this)[0] = v[0];
        (*this)[1] = v[1];
        (*this)[2] = v[2];
        (*this)[3] = v[3];
    }

    /**
     * Convert to GPU vector
     */
    operator gpu::blas::vector<T, 4>()
    {
        gpu::blas::vector<T, 4> v;
        v[0] = (*this)[0];
        v[1] = (*this)[1];
        v[2] = (*this)[2];
        v[3] = (*this)[3];
        return v;
    }
    /**
     * Convert from CUDA vector
     */
    vector(float4 const& v)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
        (*this)[2] = v.z;
        (*this)[3] = v.w;
    }

    /**
     * Convert to CUDA vector
     */
    operator float4() const
    {
        float4 v;
        v.x = (*this)[0];
        v.y = (*this)[1];
        v.z = (*this)[2];
        v.w = (*this)[3];
        return v;
    }

#endif /* WITH_CUDA */
};

template <>
struct vector<float, 4>
  : public boost::array<float, 4>
{
public:
    typedef boost::array<float, 4> _Base;
    typedef _Base::value_type value_type;
    enum { static_size = _Base::static_size };

    vector()
    {}

    /**
     * Explicit conversion from vector of convertible element type
     */
    template <typename T_>
    explicit vector(vector<T_, 4> const& v,
      typename boost::enable_if<boost::is_convertible<T_, float> >::type* dummy = 0)
    {
        (*this)[0] = static_cast<float>(v[0]);
        (*this)[1] = static_cast<float>(v[1]);
        (*this)[2] = static_cast<float>(v[2]);
        (*this)[3] = static_cast<float>(v[3]);
    }

    /**
     * Initialization by scalar
     */
    vector(float const& s)
    {
        (*this)[0] = s;
        (*this)[1] = s;
        (*this)[2] = s;
        (*this)[3] = s;
    }

    /**
     * Initialization by scalar elements
     */
    vector(float const& x, float const& y, float const& z)
    {
        (*this)[0] = x;
        (*this)[1] = y;
        (*this)[2] = z;
        (*this)[3] = z;
    }

#ifdef WITH_CUDA

    /**
     * Convert from GPU vector
     */
    template <typename T_>
    vector(gpu::blas::vector<T_, 4> const& v,
      typename boost::enable_if<boost::is_convertible<T_, float> >::type* dummy = 0)
    {
        (*this)[0] = v[0];
        (*this)[1] = v[1];
        (*this)[2] = v[2];
        (*this)[3] = v[3];
    }

    /**
     * Convert to GPU vector
     */
    operator gpu::blas::vector<float, 4>()
    {
        gpu::blas::vector<float, 4> v;
        v[0] = (*this)[0];
        v[1] = (*this)[1];
        v[2] = (*this)[2];
        v[3] = (*this)[3];
        return v;
    }
    /**
     * Convert from CUDA vector
     */
    vector(float4 const& v)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
        (*this)[2] = v.z;
        (*this)[3] = v.w;
    }

    /**
     * Convert to CUDA vector
     */
    operator float4() const
    {
        float4 v;
        v.x = (*this)[0];
        v.y = (*this)[1];
        v.z = (*this)[2];
        v.w = (*this)[3];
        return v;
    }

#endif /* WITH_CUDA */
};

/**
 * Assignment by elementwise vector addition
 */
template <typename T, typename T_>
inline typename boost::enable_if<boost::is_convertible<T_, T>, vector<T, 4>&>::type
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
inline typename boost::enable_if<boost::is_convertible<T_, T>, vector<T, 4>&>::type
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
inline typename boost::enable_if<boost::is_convertible<T_, T>, vector<T, 4>&>::type
operator*=(vector<T, 4>& v, T_ const& s)
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
inline typename boost::enable_if<boost::is_convertible<T_, T>, vector<T, 4>&>::type
operator/=(vector<T, 4>& v, T_ const& s)
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
inline typename boost::enable_if<boost::mpl::and_<boost::is_integral<T>, boost::is_integral<T_> >, vector<T, 4>&>::type
operator%=(vector<T, 4>& v, T_ const& s)
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
inline vector<T, 4> operator+(vector<T, 4> v, vector<T, 4> const& w)
{
    v[0] += w[0];
    v[1] += w[1];
    v[2] += w[2];
    v[3] += w[3];
    return v;
}

/**
 * Elementwise vector subtraction
 */
template <typename T>
inline vector<T, 4> operator-(vector<T, 4> v, vector<T, 4> const& w)
{
    v[0] -= w[0];
    v[1] -= w[1];
    v[2] -= w[2];
    v[3] -= w[3];
    return v;
}

/**
 * Scalar multiplication
 */
template <typename T, typename T_>
inline typename boost::enable_if<boost::is_convertible<T_, T>, vector<T, 4> >::type
operator*(vector<T, 4> v, T_ const& s)
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
inline typename boost::enable_if<boost::is_convertible<T_, T>, vector<T, 4> >::type
operator*(T_ const& s, vector<T, 4> v)
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
inline typename boost::enable_if<boost::is_convertible<T_, T>, vector<T, 4> >::type
operator/(vector<T, 4> v, T_ const& s)
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
inline typename boost::enable_if<boost::mpl::and_<boost::is_integral<T>, boost::is_integral<T_> >, vector<T, 4> >::type
operator%(vector<T, 4> v, T_ const& s)
{
    v[0] %= s;
    v[1] %= s;
    v[2] %= s;
    v[3] %= s;
    return v;
}

/**
 * Unary negative operator
 */
template <typename T>
inline vector<T, 4> operator-(vector<T, 4> v)
{
    v[0] = -v[0];
    v[1] = -v[1];
    v[2] = -v[2];
    v[3] = -v[3];
    return v;
}

/**
 * Inner product
 */
template <typename T>
inline T inner_prod(vector<T, 4> const& v, vector<T, 4> const& w)
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
inline vector<T, 4> element_prod(vector<T, 4> v, vector<T, 4> const& w)
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
inline vector<T, 4> element_div(vector<T, 4> v, vector<T, 4> const& w)
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
inline typename boost::enable_if<boost::is_integral<T>, vector<T, 4> >::type
element_mod(vector<T, 4> v, vector<T, 4> const& w)
{
    v[0] %= w[0];
    v[1] %= w[1];
    v[2] %= w[2];
    v[3] %= w[3];
    return v;
}

/**
 * Write vector elements to output stream
 */
template <typename T>
inline std::ostream& operator<<(std::ostream& os, vector<T, 4> const& v)
{
    os << v[0] << " " << v[1] << " " << v[2] << " " << v[3];
    return os;
}

/**
 * Read vector elements from input stream
 */
template <typename T>
inline std::istream& operator>>(std::istream& is, vector<T, 4>& v)
{
    is >> v[0] >> v[1] >> v[2] >> v[3];
    return is;
}

/**
 * Elementwise round to nearest integer
 */
inline vector<float, 4> rint(vector<float, 4> v)
{
    v[0] = ::rintf(v[0]);
    v[1] = ::rintf(v[1]);
    v[2] = ::rintf(v[2]);
    v[3] = ::rintf(v[3]);
    return v;
}

inline vector<double, 4> rint(vector<double, 4> v)
{
    v[0] = ::rint(v[0]);
    v[1] = ::rint(v[1]);
    v[2] = ::rint(v[2]);
    v[3] = ::rint(v[3]);
    return v;
}

/**
 * Elementwise round to nearest integer, away from zero
 */
inline vector<float, 4> round(vector<float, 4> v)
{
    v[0] = ::roundf(v[0]);
    v[1] = ::roundf(v[1]);
    v[2] = ::roundf(v[2]);
    v[3] = ::roundf(v[3]);
    return v;
}

inline vector<double, 4> round(vector<double, 4> v)
{
    v[0] = ::round(v[0]);
    v[1] = ::round(v[1]);
    v[2] = ::round(v[2]);
    v[3] = ::round(v[3]);
    return v;
}

/**
 * Elementwise round to nearest integer not greater than argument
 */
inline vector<float, 4> floor(vector<float, 4> v)
{
    v[0] = std::floor(v[0]);
    v[1] = std::floor(v[1]);
    v[2] = std::floor(v[2]);
    v[3] = std::floor(v[3]);
    return v;
}

inline vector<double, 4> floor(vector<double, 4> v)
{
    v[0] = std::floor(v[0]);
    v[1] = std::floor(v[1]);
    v[2] = std::floor(v[2]);
    v[3] = std::floor(v[3]);
    return v;
}

/**
 * Elementwise round to nearest integer not less argument
 */
inline vector<float, 4> ceil(vector<float, 4> v)
{
    v[0] = std::ceil(v[0]);
    v[1] = std::ceil(v[1]);
    v[2] = std::ceil(v[2]);
    v[3] = std::ceil(v[3]);
    return v;
}

inline vector<double, 4> ceil(vector<double, 4> v)
{
    v[0] = std::ceil(v[0]);
    v[1] = std::ceil(v[1]);
    v[2] = std::ceil(v[2]);
    v[3] = std::ceil(v[3]);
    return v;
}

/**
 * Elementwise square root function
 */
inline vector<float, 4> sqrt(vector<float, 4> v)
{
    v[0] = std::sqrt(v[0]);
    v[1] = std::sqrt(v[1]);
    v[2] = std::sqrt(v[2]);
    v[3] = std::sqrt(v[3]);
    return v;
}

inline vector<double, 4> sqrt(vector<double, 4> v)
{
    v[0] = std::sqrt(v[0]);
    v[1] = std::sqrt(v[1]);
    v[2] = std::sqrt(v[2]);
    v[3] = std::sqrt(v[3]);
    return v;
}

/**
 * Elementwise cos function
 */
inline vector<float, 4> cos(vector<float, 4> v)
{
    v[0] = std::cos(v[0]);
    v[1] = std::cos(v[1]);
    v[2] = std::cos(v[2]);
    v[3] = std::cos(v[3]);
    return v;
}

inline vector<double, 4> cos(vector<double, 4> v)
{
    v[0] = std::cos(v[0]);
    v[1] = std::cos(v[1]);
    v[2] = std::cos(v[2]);
    v[3] = std::cos(v[3]);
    return v;
}

/**
 * Elementwise sin function
 */
inline vector<float, 4> sin(vector<float, 4> v)
{
    v[0] = std::sin(v[0]);
    v[1] = std::sin(v[1]);
    v[2] = std::sin(v[2]);
    v[3] = std::sin(v[3]);
    return v;
}

inline vector<double, 4> sin(vector<double, 4> v)
{
    v[0] = std::sin(v[0]);
    v[1] = std::sin(v[1]);
    v[2] = std::sin(v[2]);
    v[3] = std::sin(v[3]);
    return v;
}

} // namespace detail

}}}} // namespace halmd::numeric::host::blas

#endif /* ! HALMD_NUMERIC_HOST_BLAS_DETAIL_VECTOR4D_HPP */
