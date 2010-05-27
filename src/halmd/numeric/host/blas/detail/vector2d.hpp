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

#ifndef HALMD_NUMERIC_HOST_BLAS_DETAIL_VECTOR2D_HPP
#define HALMD_NUMERIC_HOST_BLAS_DETAIL_VECTOR2D_HPP

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
 * Two-dimensional floating-point vector
 */
template <typename T>
struct vector<T, 2>
  : public boost::array<T, 2>
{
public:
    typedef boost::array<T, 2> _Base;
    typedef typename _Base::value_type value_type;
    enum { static_size = _Base::static_size };

    vector()
    {}

    /**
     * Explicit conversion from vector of convertible element type
     */
    template <typename T_>
    explicit vector(vector<T_, 2> const& v,
      typename boost::enable_if<boost::is_convertible<T_, T> >::type* dummy = 0)
    {
        (*this)[0] = static_cast<T>(v[0]);
        (*this)[1] = static_cast<T>(v[1]);
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
    }

#ifdef WITH_CUDA

    /**
     * Convert from GPU vector
     */
    template <typename T_>
    vector(gpu::blas::vector<T_, 2> const& v,
      typename boost::enable_if<boost::is_convertible<T_, T> >::type* dummy = 0)
    {
        (*this)[0] = v[0];
        (*this)[1] = v[1];
    }

    /**
     * Convert to GPU vector
     */
    operator gpu::blas::vector<T, 2>()
    {
        gpu::blas::vector<T, 2> v;
        v[0] = (*this)[0];
        v[1] = (*this)[1];
        return v;
    }

#endif /* WITH_CUDA */
};

template <>
struct vector<float, 2>
  : public boost::array<float, 2>
{
public:
    typedef boost::array<float, 2> _Base;
    typedef _Base::value_type value_type;
    enum { static_size = _Base::static_size };

    vector()
    {}

    /**
     * Explicit conversion from vector of convertible element type
     */
    template <typename T_>
    explicit vector(vector<T_, 2> const& v,
      typename boost::enable_if<boost::is_convertible<T_, float> >::type* dummy = 0)
    {
        (*this)[0] = static_cast<float>(v[0]);
        (*this)[1] = static_cast<float>(v[1]);
    }

    /**
     * Initialization by scalar
     */
    vector(float const& s)
    {
        (*this)[0] = s;
        (*this)[1] = s;
    }

    /**
     * Initialization by scalar elements
     */
    vector(float const& x, float const& y, float const& z)
    {
        (*this)[0] = x;
        (*this)[1] = y;
    }

#ifdef WITH_CUDA

    /**
     * Convert from GPU vector
     */
    template <typename T_>
    vector(gpu::blas::vector<T_, 2> const& v,
      typename boost::enable_if<boost::is_convertible<T_, float> >::type* dummy = 0)
    {
        (*this)[0] = v[0];
        (*this)[1] = v[1];
    }

    /**
     * Convert to GPU vector
     */
    operator gpu::blas::vector<float, 2>()
    {
        gpu::blas::vector<float, 2> v;
        v[0] = (*this)[0];
        v[1] = (*this)[1];
        return v;
    }

    /**
     * Convert from CUDA vector
     */
    vector(float2 const& v)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
    }

    vector(float3 const& v)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
    }

    vector(float4 const& v)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
    }

    /**
     * Convert to CUDA vector
     */
    operator float2() const
    {
        float2 v;
        v.x = (*this)[0];
        v.y = (*this)[1];
        return v;
    }

    operator float3() const
    {
        float3 v;
        v.x = (*this)[0];
        v.y = (*this)[1];
        return v;
    }

    operator float4() const
    {
        float4 v;
        v.x = (*this)[0];
        v.y = (*this)[1];
        return v;
    }

#endif /* WITH_CUDA */
};

/**
 * Assignment by elementwise vector addition
 */
template <typename T, typename T_>
inline typename boost::enable_if<boost::is_convertible<T_, T>, vector<T, 2>&>::type
operator+=(vector<T, 2>& v, vector<T_, 2> const& w)
{
    v[0] += w[0];
    v[1] += w[1];
    return v;
}

/**
 * Assignment by elementwise vector subtraction
 */
template <typename T, typename T_>
inline typename boost::enable_if<boost::is_convertible<T_, T>, vector<T, 2>&>::type
operator-=(vector<T, 2>& v, vector<T_, 2> const& w)
{
    v[0] -= w[0];
    v[1] -= w[1];
    return v;
}

/**
 * Assignment by scalar multiplication
 */
template <typename T, typename T_>
inline typename boost::enable_if<boost::is_convertible<T_, T>, vector<T, 2>&>::type
operator*=(vector<T, 2>& v, T_ const& s)
{
    v[0] *= s;
    v[1] *= s;
    return v;
}

/**
 * Assignment by scalar division
 */
template <typename T, typename T_>
inline typename boost::enable_if<boost::is_convertible<T_, T>, vector<T, 2>&>::type
operator/=(vector<T, 2>& v, T_ const& s)
{
    v[0] /= s;
    v[1] /= s;
    return v;
}

/**
 * Assignment by scalar modulus
 */
template <typename T, typename T_>
inline typename boost::enable_if<boost::mpl::and_<boost::is_integral<T>, boost::is_integral<T_> >, vector<T, 2>&>::type
operator%=(vector<T, 2>& v, T_ const& s)
{
    v[0] %= s;
    v[1] %= s;
    return v;
}

/**
 * Elementwise vector addition
 */
template <typename T>
inline vector<T, 2> operator+(vector<T, 2> v, vector<T, 2> const& w)
{
    v[0] += w[0];
    v[1] += w[1];
    return v;
}

/**
 * Elementwise vector subtraction
 */
template <typename T>
inline vector<T, 2> operator-(vector<T, 2> v, vector<T, 2> const& w)
{
    v[0] -= w[0];
    v[1] -= w[1];
    return v;
}

/**
 * Scalar multiplication
 */
template <typename T, typename T_>
inline typename boost::enable_if<boost::is_convertible<T_, T>, vector<T, 2> >::type
operator*(vector<T, 2> v, T_ const& s)
{
    v[0] *= s;
    v[1] *= s;
    return v;
}

/**
 * Scalar multiplication
 */
template <typename T, typename T_>
inline typename boost::enable_if<boost::is_convertible<T_, T>, vector<T, 2> >::type
operator*(T_ const& s, vector<T, 2> v)
{
    v[0] *= s;
    v[1] *= s;
    return v;
}

/**
 * Scalar division
 */
template <typename T, typename T_>
inline typename boost::enable_if<boost::is_convertible<T_, T>, vector<T, 2> >::type
operator/(vector<T, 2> v, T_ const& s)
{
    v[0] /= s;
    v[1] /= s;
    return v;
}

/**
 * Scalar modulus
 */
template <typename T, typename T_>
inline typename boost::enable_if<boost::mpl::and_<boost::is_integral<T>, boost::is_integral<T_> >, vector<T, 2> >::type
operator%(vector<T, 2> v, T_ const& s)
{
    v[0] %= s;
    v[1] %= s;
    return v;
}

/**
 * Unary negative operator
 */
template <typename T>
inline vector<T, 2> operator-(vector<T, 2> v)
{
    v[0] = -v[0];
    v[1] = -v[1];
    return v;
}

/**
 * Inner product
 */
template <typename T>
inline T inner_prod(vector<T, 2> const& v, vector<T, 2> const& w)
{
    T s = v[0] * w[0];
    s  += v[1] * w[1];
    return s;
}

/**
 * Elementwise vector multiplication
 */
template <typename T>
inline vector<T, 2> element_prod(vector<T, 2> v, vector<T, 2> const& w)
{
    v[0] *= w[0];
    v[1] *= w[1];
    return v;
}

/**
 * Elementwise vector division
 */
template <typename T>
inline vector<T, 2> element_div(vector<T, 2> v, vector<T, 2> const& w)
{
    v[0] /= w[0];
    v[1] /= w[1];
    return v;
}

/**
 * Elementwise vector modulus
 */
template <typename T>
inline typename boost::enable_if<boost::is_integral<T>, vector<T, 2> >::type
element_mod(vector<T, 2> v, vector<T, 2> const& w)
{
    v[0] %= w[0];
    v[1] %= w[1];
    return v;
}

/**
 * Write vector elements to output stream
 */
template <typename T>
inline std::ostream& operator<<(std::ostream& os, vector<T, 2> const& v)
{
    os << v[0] << " " << v[1];
    return os;
}

/**
 * Read vector elements from input stream
 */
template <typename T>
inline std::istream& operator>>(std::istream& is, vector<T, 2>& v)
{
    is >> v[0] >> v[1];
    return is;
}

/**
 * Elementwise round to nearest integer
 */
inline vector<float, 2> rint(vector<float, 2> v)
{
    v[0] = ::rintf(v[0]);
    v[1] = ::rintf(v[1]);
    return v;
}

inline vector<double, 2> rint(vector<double, 2> v)
{
    v[0] = ::rint(v[0]);
    v[1] = ::rint(v[1]);
    return v;
}

/**
 * Elementwise round to nearest integer, away from zero
 */
inline vector<float, 2> round(vector<float, 2> v)
{
    v[0] = ::roundf(v[0]);
    v[1] = ::roundf(v[1]);
    return v;
}

inline vector<double, 2> round(vector<double, 2> v)
{
    v[0] = ::round(v[0]);
    v[1] = ::round(v[1]);
    return v;
}

/**
 * Elementwise round to nearest integer not greater than argument
 */
inline vector<float, 2> floor(vector<float, 2> v)
{
    v[0] = std::floor(v[0]);
    v[1] = std::floor(v[1]);
    return v;
}

inline vector<double, 2> floor(vector<double, 2> v)
{
    v[0] = std::floor(v[0]);
    v[1] = std::floor(v[1]);
    return v;
}

/**
 * Elementwise round to nearest integer not less argument
 */
inline vector<float, 2> ceil(vector<float, 2> v)
{
    v[0] = std::ceil(v[0]);
    v[1] = std::ceil(v[1]);
    return v;
}

inline vector<double, 2> ceil(vector<double, 2> v)
{
    v[0] = std::ceil(v[0]);
    v[1] = std::ceil(v[1]);
    return v;
}

/**
 * Elementwise square root function
 */
inline vector<float, 2> sqrt(vector<float, 2> v)
{
    v[0] = std::sqrt(v[0]);
    v[1] = std::sqrt(v[1]);
    return v;
}

inline vector<double, 2> sqrt(vector<double, 2> v)
{
    v[0] = std::sqrt(v[0]);
    v[1] = std::sqrt(v[1]);
    return v;
}

/**
 * Elementwise cos function
 */
inline vector<float, 2> cos(vector<float, 2> v)
{
    v[0] = std::cos(v[0]);
    v[1] = std::cos(v[1]);
    return v;
}

inline vector<double, 2> cos(vector<double, 2> v)
{
    v[0] = std::cos(v[0]);
    v[1] = std::cos(v[1]);
    return v;
}

/**
 * Elementwise sin function
 */
inline vector<float, 2> sin(vector<float, 2> v)
{
    v[0] = std::sin(v[0]);
    v[1] = std::sin(v[1]);
    return v;
}

inline vector<double, 2> sin(vector<double, 2> v)
{
    v[0] = std::sin(v[0]);
    v[1] = std::sin(v[1]);
    return v;
}

} // namespace detail

}}}} // namespace halmd::numeric::host::blas

#endif /* ! HALMD_NUMERIC_HOST_BLAS_DETAIL_VECTOR2D_HPP */
