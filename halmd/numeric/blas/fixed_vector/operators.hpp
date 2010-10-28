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

#ifndef HALMD_NUMERIC_BLAS_FIXED_VECTOR_OPERATORS_HPP
#define HALMD_NUMERIC_BLAS_FIXED_VECTOR_OPERATORS_HPP

#include <boost/mpl/and.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/type_traits/is_floating_point.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/utility/enable_if.hpp>
#ifndef __CUDACC__
# include <cmath>
# include <iostream>
#endif
#ifdef WITH_CUDA
# include <cuda_runtime.h> // CUDA vector types for host compiler
#endif

#include <halmd/config.hpp>

namespace halmd
{
namespace detail { namespace numeric { namespace blas
{

HALMD_GPU_USING(algorithm::gpu::tuple, boost::tuple);
HALMD_GPU_USING(algorithm::gpu::tie, boost::tie);
HALMD_GPU_USING(algorithm::gpu::make_tuple, boost::make_tuple);

/**
 * Returns "high" and "low" single precision vector tuple
 */
template <size_t N>
inline HALMD_GPU_ENABLED
tuple<fixed_vector<float, N>, fixed_vector<float, N> >
split(fixed_vector<dsfloat, N> const& v)
{
    fixed_vector<float, N> hi, lo;
    for (size_t i = 0; i < N; ++i) {
        tie(hi[i], lo[i]) = split(v[i]);
    }
    return make_tuple(hi, lo);
}

/**
 * Assignment by elementwise vector addition
 */
template <typename T, typename U, size_t N>
inline HALMD_GPU_ENABLED
typename boost::enable_if<boost::is_convertible<U, T>, fixed_vector<T, N>&>::type
operator+=(fixed_vector<T, N>& v, fixed_vector<U, N> const& w)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] += w[i];
    }
    return v;
}

/**
 * Assignment by elementwise vector subtraction
 */
template <typename T, typename U, size_t N>
inline HALMD_GPU_ENABLED
typename boost::enable_if<boost::is_convertible<U, T>, fixed_vector<T, N>&>::type
operator-=(fixed_vector<T, N>& v, fixed_vector<U, N> const& w)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] -= w[i];
    }
    return v;
}

/**
 * Assignment by scalar multiplication
 */
template <typename T, typename U, size_t N>
inline HALMD_GPU_ENABLED
typename boost::enable_if<boost::is_convertible<U, T>, fixed_vector<T, N>&>::type
operator*=(fixed_vector<T, N>& v, U s)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] *= s;
    }
    return v;
}

/**
 * Assignment by scalar division
 */
template <typename T, typename U, size_t N>
inline HALMD_GPU_ENABLED
typename boost::enable_if<boost::is_convertible<U, T>, fixed_vector<T, N>&>::type
operator/=(fixed_vector<T, N>& v, U s)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] /= s;
    }
    return v;
}

/**
 * Assignment by scalar modulus
 */
template <typename T, typename U, size_t N>
inline HALMD_GPU_ENABLED
typename boost::enable_if<boost::mpl::and_<boost::is_integral<T>, boost::is_integral<U> >, fixed_vector<T, N>&>::type
operator%=(fixed_vector<T, N>& v, U s)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] %= s;
    }
    return v;
}

/**
 * Elementwise vector addition
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<T, N> operator+(fixed_vector<T, N> v, fixed_vector<T, N> const& w)
{
    v += w;
    return v;
}

/**
 * Elementwise vector subtraction
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<T, N> operator-(fixed_vector<T, N> v, fixed_vector<T, N> const& w)
{
    v -= w;
    return v;
}

/**
 * Elementwise change of sign
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<T, N> operator-(fixed_vector<T, N> v)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] = -v[i];
    }
    return v;
}

/**
 * Scalar multiplication
 */
template <typename T, typename U, size_t N>
inline HALMD_GPU_ENABLED
typename boost::enable_if<boost::is_convertible<U, T>, fixed_vector<T, N> >::type
operator*(fixed_vector<T, N> v, U s)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] *= s;
    }
    return v;
}

/**
 * Scalar multiplication
 */
template <typename T, typename U, size_t N>
inline HALMD_GPU_ENABLED
typename boost::enable_if<boost::is_convertible<U, T>, fixed_vector<T, N> >::type
operator*(U s, fixed_vector<T, N> v)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] *= s;
    }
    return v;
}

/**
 * Scalar division
 */
template <typename T, typename U, size_t N>
inline HALMD_GPU_ENABLED
typename boost::enable_if<boost::is_convertible<U, T>, fixed_vector<T, N> >::type
operator/(fixed_vector<T, N> v, U s)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] /= s;
    }
    return v;
}

/**
 * Scalar modulus
 */
template <typename T, typename U, size_t N>
inline HALMD_GPU_ENABLED
typename boost::enable_if<boost::mpl::and_<boost::is_integral<T>, boost::is_integral<U> >, fixed_vector<T, N> >::type
operator%(fixed_vector<T, N> v, U s)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] %= s;
    }
    return v;
}

/**
 * Inner product
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED
T inner_prod(fixed_vector<T, N> const& v, fixed_vector<T, N> const& w)
{
    T s = v[0] * w[0];
    for (size_t i = 1; i < N; ++i) {
        s += v[i] * w[i];
    }
    return s;
}

/**
 * Elementwise vector multiplication
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<T, N> element_prod(fixed_vector<T, N> v, fixed_vector<T, N> const& w)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] *= w[i];
    }
    return v;
}

/**
 * Elementwise vector division
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<T, N> element_div(fixed_vector<T, N> v, fixed_vector<T, N> const& w)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] /= w[i];
    }
    return v;
}

/**
 * Elementwise vector modulus
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED
typename boost::enable_if<boost::is_integral<T>, fixed_vector<T, N> >::type
element_mod(fixed_vector<T, N> v, fixed_vector<T, N> const& w)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] %= w[i];
    }
    return v;
}

/**
 * Elementwise round to nearest integer not greater than argument
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED
typename boost::enable_if<boost::is_floating_point<T>, fixed_vector<T, N> >::type
floor(fixed_vector<T, N> v)
{
    HALMD_GPU_USING(::floor, std::floor);
    for (size_t i = 0; i < N; ++i) {
        v[i] = floor(v[i]);
    }
    return v;
}

/**
 * Elementwise round to nearest integer not less argument
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED
typename boost::enable_if<boost::is_floating_point<T>, fixed_vector<T, N> >::type
ceil(fixed_vector<T, N> v)
{
    HALMD_GPU_USING(::ceil, std::ceil);
    for (size_t i = 0; i < N; ++i) {
        v[i] = ceil(v[i]);
    }
    return v;
}

/**
 * Elementwise square root function
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED
typename boost::enable_if<boost::is_floating_point<T>, fixed_vector<T, N> >::type
sqrt(fixed_vector<T, N> v)
{
    HALMD_GPU_USING(::sqrt, std::sqrt);
    for (size_t i = 0; i < N; ++i) {
        v[i] = sqrt(v[i]);
    }
    return v;
}

/**
 * Elementwise cosine function
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED
typename boost::enable_if<boost::is_floating_point<T>, fixed_vector<T, N> >::type
cos(fixed_vector<T, N> v)
{
    HALMD_GPU_USING(::cos, std::cos);
    for (size_t i = 0; i < N; ++i) {
        v[i] = cos(v[i]);
    }
    return v;
}

/**
 * Elementwise sine function
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED
typename boost::enable_if<boost::is_floating_point<T>, fixed_vector<T, N> >::type
sin(fixed_vector<T, N> v)
{
    HALMD_GPU_USING(::sin, std::sin);
    for (size_t i = 0; i < N; ++i) {
        v[i] = sin(v[i]);
    }
    return v;
}

/**
 * Elementwise absolute value
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED
typename boost::enable_if<boost::is_floating_point<T>, fixed_vector<T, N> >::type
fabs(fixed_vector<T, N> v)
{
    HALMD_GPU_USING(::fabs, std::fabs);
    for (size_t i = 0; i < N; ++i) {
        v[i] = fabs(v[i]);
    }
    return v;
}

/**
 * Floating-point remainder function, round towards nearest integer
 */
template <size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<double, N>
remainder(fixed_vector<double, N> v, fixed_vector<double, N> const& w)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] = ::remainder(v[i], w[i]);
    }
    return v;
}

template <size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<float, N>
remainder(fixed_vector<float, N> v, fixed_vector<float, N> const& w)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] = ::remainderf(v[i], w[i]);
    }
    return v;
}

/**
 * Floating-point remainder function, round towards zero
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED
typename boost::enable_if<boost::is_floating_point<T>, fixed_vector<T, N> >::type
fmod(fixed_vector<T, N> v, fixed_vector<T, N> const& w)
{
    HALMD_GPU_USING(::fmod, std::fmod);
    for (size_t i = 0; i < N; ++i) {
        v[i] = fmod(v[i], w[i]);
    }
    return v;
}

/**
 * Elementwise round to nearest integer
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED
typename boost::enable_if<boost::is_same<T, double>, fixed_vector<T, N> >::type
rint(fixed_vector<T, N> v)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] = ::rint(v[i]);
    }
    return v;
}

template <typename T, size_t N>
inline HALMD_GPU_ENABLED
typename boost::enable_if<boost::is_same<T, float>, fixed_vector<T, N> >::type
rint(fixed_vector<T, N> v)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] = ::rintf(v[i]);
    }
    return v;
}

/**
 * Elementwise round to nearest integer, away from zero
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED
typename boost::enable_if<boost::is_same<T, double>, fixed_vector<T, N> >::type
round(fixed_vector<T, N> v)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] = ::round(v[i]);
    }
    return v;
}

template <typename T, size_t N>
inline HALMD_GPU_ENABLED
typename boost::enable_if<boost::is_same<T, float>, fixed_vector<T, N> >::type
round(fixed_vector<T, N> v)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] = ::roundf(v[i]);
    }
    return v;
}

#ifndef __CUDACC__

/**
 * Write vector elements to output stream
 */
template <typename T, size_t N>
inline std::ostream& operator<<(std::ostream& os, fixed_vector<T, N> const& v)
{
    os << v[0];
    for (size_t i = 1; i < N; ++i) {
        os << " " << v[i];
    }
    return os;
}

/**
 * Read vector elements from input stream
 */
template <typename T, size_t N>
inline std::istream& operator>>(std::istream& is, fixed_vector<T, N>& v)
{
    for (size_t i = 0; i < N; ++i) {
        is >> v[i];
    }
    return is;
}

#endif /* ! __CUDACC__ */

}}} // namespace detail::numeric::blas

// import into top-level namespace
using detail::numeric::blas::fixed_vector;

} // namespace halmd

#endif /* ! HALMD_NUMERIC_BLAS_FIXED_VECTOR_OPERATORS_HPP */
