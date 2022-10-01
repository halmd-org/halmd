/*
 * Copyright © 2008-2013 Felix Höfling
 * Copyright © 2008-2012 Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#ifndef HALMD_NUMERIC_BLAS_DETAIL_OPERATORS_HPP
#define HALMD_NUMERIC_BLAS_DETAIL_OPERATORS_HPP

#include <halmd/config.hpp>

#include <algorithm>
#include <type_traits>

#ifndef __CUDACC__
# include <cmath>
# include <iostream>
#endif

// CUDA vector types for host compiler
#ifdef HALMD_WITH_GPU
# include <cuda_runtime.h>
#endif

#include <halmd/numeric/blas/detail/vector.hpp>
#include <halmd/utility/tuple.hpp>

namespace halmd {
namespace numeric {
namespace blas {
namespace detail {

/**
 * Returns "high" and "low" single precision vector tuple
 */
template <size_t N, size_t L, size_t U>
inline HALMD_GPU_ENABLED
typename std::enable_if<(U <= L), void>::type
split(fixed_vector<dsfloat, N> const& v, fixed_vector<float, N>& hi, fixed_vector<float, N>& lo)
{
    tie(hi[L], lo[L]) = split(v[L]);
}

template <size_t N, size_t L, size_t U>
inline HALMD_GPU_ENABLED
typename std::enable_if<(U > L), void>::type
split(fixed_vector<dsfloat, N> const& v, fixed_vector<float, N>& hi, fixed_vector<float, N>& lo)
{
    split<N, L, (L + U) / 2>(v, hi, lo);
    split<N, (L + U) / 2 + 1, U>(v, hi, lo);
}

template <size_t N>
inline HALMD_GPU_ENABLED
tuple<fixed_vector<float, N>, fixed_vector<float, N> >
split(fixed_vector<dsfloat, N> const& v)
{
    fixed_vector<float, N> hi, lo;
    split<N, 0, N - 1>(v, hi, lo);
    return make_tuple(hi, lo);
}

/**
 * Assignment by elementwise vector addition
 */
template <typename T, typename S, size_t N, size_t L, size_t U>
inline HALMD_GPU_ENABLED
typename std::enable_if<(U <= L), void>::type
operator+=(fixed_vector<T, N>& v, fixed_vector<S, N> const& w)
{
    v[L] += static_cast<T>(w[L]);
}

template <typename T, typename S, size_t N, size_t L, size_t U>
inline HALMD_GPU_ENABLED
typename std::enable_if<(U > L), void>::type
operator+=(fixed_vector<T, N>& v, fixed_vector<S, N> const& w)
{
    operator+=<T, S, N, L, (L + U) / 2>(v, w);
    operator+=<T, S, N, (L + U) / 2 + 1, U>(v, w);
}

template <typename T, typename S, size_t N>
inline HALMD_GPU_ENABLED
typename std::enable_if<std::is_convertible<S, T>::value, fixed_vector<T, N>&>::type
operator+=(fixed_vector<T, N>& v, fixed_vector<S, N> const& w)
{
    operator+=<T, S, N, 0, N - 1>(v, w);
    return v;
}

/**
 * Assignment by elementwise vector subtraction
 */
template <typename T, typename S, size_t N, size_t L, size_t U>
inline HALMD_GPU_ENABLED
typename std::enable_if<(U <= L), void>::type
operator-=(fixed_vector<T, N>& v, fixed_vector<S, N> const& w)
{
    v[L] -= static_cast<T>(w[L]);
}

template <typename T, typename S, size_t N, size_t L, size_t U>
inline HALMD_GPU_ENABLED
typename std::enable_if<(U > L), void>::type
operator-=(fixed_vector<T, N>& v, fixed_vector<S, N> const& w)
{
    operator-=<T, S, N, L, (L + U) / 2>(v, w);
    operator-=<T, S, N, (L + U) / 2 + 1, U>(v, w);
}

template <typename T, typename S, size_t N>
inline HALMD_GPU_ENABLED
typename std::enable_if<std::is_convertible<S, T>::value, fixed_vector<T, N>&>::type
operator-=(fixed_vector<T, N>& v, fixed_vector<S, N> const& w)
{
    operator-=<T, S, N, 0, N - 1>(v, w);
    return v;
}

/**
 * Assignment by scalar multiplication
 */
template <typename T, typename S, size_t N, size_t L, size_t U>
inline HALMD_GPU_ENABLED
typename std::enable_if<(U <= L), void>::type
operator*=(fixed_vector<T, N>& v, S const& s)
{
    v[L] *= static_cast<T>(s);
}

template <typename T, typename S, size_t N, size_t L, size_t U>
inline HALMD_GPU_ENABLED
typename std::enable_if<(U > L), void>::type
operator*=(fixed_vector<T, N>& v, S const& s)
{
    operator*=<T, S, N, L, (L + U) / 2>(v, s);
    operator*=<T, S, N, (L + U) / 2 + 1, U>(v, s);
}

template <typename T, typename S, size_t N>
inline HALMD_GPU_ENABLED
typename std::enable_if<std::is_convertible<S, T>::value, fixed_vector<T, N>&>::type
operator*=(fixed_vector<T, N>& v, S const& s)
{
    operator*=<T, S, N, 0, N - 1>(v, s);
    return v;
}

/**
 * Assignment by scalar division
 */
template <typename T, typename S, size_t N, size_t L, size_t U>
inline HALMD_GPU_ENABLED
typename std::enable_if<(U <= L), void>::type
operator/=(fixed_vector<T, N>& v, S const& s)
{
    v[L] /= static_cast<T>(s);
}

template <typename T, typename S, size_t N, size_t L, size_t U>
inline HALMD_GPU_ENABLED
typename std::enable_if<(U > L), void>::type
operator/=(fixed_vector<T, N>& v, S const& s)
{
    operator/=<T, S, N, L, (L + U) / 2>(v, s);
    operator/=<T, S, N, (L + U) / 2 + 1, U>(v, s);
}

template <typename T, typename S, size_t N>
inline HALMD_GPU_ENABLED
typename std::enable_if<std::is_convertible<S, T>::value, fixed_vector<T, N>&>::type
operator/=(fixed_vector<T, N>& v, S const& s)
{
    operator/=<T, S, N, 0, N - 1>(v, s);
    return v;
}

/**
 * Assignment by scalar modulus
 */
template <typename T, typename S, size_t N, size_t L, size_t U>
inline HALMD_GPU_ENABLED
typename std::enable_if<(U <= L), void>::type
operator%=(fixed_vector<T, N>& v, S const& s)
{
    v[L] %= static_cast<T>(s);
}

template <typename T, typename S, size_t N, size_t L, size_t U>
inline HALMD_GPU_ENABLED
typename std::enable_if<(U > L), void>::type
operator%=(fixed_vector<T, N>& v, S const& s)
{
    operator%=<T, S, N, L, (L + U) / 2>(v, s);
    operator%=<T, S, N, (L + U) / 2 + 1, U>(v, s);
}

template <typename T, typename S, size_t N>
inline HALMD_GPU_ENABLED
typename std::enable_if<std::is_integral<T>::value && std::is_integral<S>::value, fixed_vector<T, N>&>::type
operator%=(fixed_vector<T, N>& v, S const& s)
{
    operator%=<T, S, N, 0, N - 1>(v, s);
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
template <typename T, size_t N, size_t L, size_t U>
inline HALMD_GPU_ENABLED
typename std::enable_if<(U <= L), void>::type
operator-(fixed_vector<T, N>& v)
{
    v[L] = -v[L];
}

template <typename T, size_t N, size_t L, size_t U>
inline HALMD_GPU_ENABLED
typename std::enable_if<(U > L), void>::type
operator-(fixed_vector<T, N>& v)
{
    operator-<T, N, L, (L + U) / 2>(v);
    operator-<T, N, (L + U) / 2 + 1, U>(v);
}

template <typename T, size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<T, N> operator-(fixed_vector<T, N> v)
{
    operator-<T, N, 0, N - 1>(v);
    return v;
}

/**
 * Scalar multiplication
 */
template <typename T, typename S, size_t N>
inline HALMD_GPU_ENABLED
typename std::enable_if<std::is_convertible<S, T>::value, fixed_vector<T, N>>::type
operator*(fixed_vector<T, N> v, S const& s)
{
    v *= s;
    return v;
}

/**
 * Scalar multiplication
 */
template <typename T, typename S, size_t N>
inline HALMD_GPU_ENABLED
typename std::enable_if<std::is_convertible<S, T>::value, fixed_vector<T, N>>::type
operator*(S const& s, fixed_vector<T, N> v)
{
    v *= s;
    return v;
}

/**
 * Scalar division
 */
template <typename T, typename S, size_t N>
inline HALMD_GPU_ENABLED
typename std::enable_if<std::is_convertible<S, T>::value, fixed_vector<T, N>>::type
operator/(fixed_vector<T, N> v, S const& s)
{
    v /= s;
    return v;
}

/**
 * Scalar modulus
 */
template <typename T, typename S, size_t N>
inline HALMD_GPU_ENABLED
typename std::enable_if<std::is_integral<T>::value && std::is_integral<S>::value, fixed_vector<T, N>>::type
operator%(fixed_vector<T, N> v, S const& s)
{
    v %= s;
    return v;
}

/**
 * Elementwise comparison
 * @returns true if _all_ elements are equal
 */
template <typename T, typename S, size_t N, size_t L, size_t U>
inline HALMD_GPU_ENABLED
typename std::enable_if<(U <= L), bool>::type
operator==(fixed_vector<T, N> const& v, fixed_vector<S, N> const& w)
{
    return v[L] == static_cast<T>(w[L]);
}

template <typename T, typename S, size_t N, size_t L, size_t U>
inline HALMD_GPU_ENABLED
typename std::enable_if<(U > L), bool>::type
operator==(fixed_vector<T, N> const& v, fixed_vector<S, N> const& w)
{
    return
        operator==<T, S, N, L, (L + U) / 2>(v, w)
     && operator==<T, S, N, (L + U) / 2 + 1, U>(v, w);
}

template <typename T, typename S, size_t N>
inline HALMD_GPU_ENABLED
typename std::enable_if<std::is_convertible<S, T>::value, bool>::type
operator==(fixed_vector<T, N> const& v, fixed_vector<S, N> const& w)
{
    return operator==<T, S, N, 0, N - 1>(v, w);
}

/**
 * Inner product
 */
template <typename T, size_t N, size_t L, size_t S>
inline HALMD_GPU_ENABLED
typename std::enable_if<(S <= L), T>::type
inner_prod(fixed_vector<T, N> const& v, fixed_vector<T, N> const& w)
{
    return v[L] * w[L];
}

template <typename T, size_t N, size_t L, size_t S>
inline HALMD_GPU_ENABLED
typename std::enable_if<(S > L), T>::type
inner_prod(fixed_vector<T, N> const& v, fixed_vector<T, N> const& w)
{
    return inner_prod<T, N, L, (L + S) / 2>(v, w) + inner_prod<T, N, (L + S) / 2 + 1, S>(v, w);
}

template <typename T, size_t N>
inline HALMD_GPU_ENABLED
T inner_prod(fixed_vector<T, N> const& v, fixed_vector<T, N> const& w)
{
    return inner_prod<T, N, 0, N - 1>(v, w);
}

/**
 * Elementwise vector multiplication
 */
template <typename T, size_t N, size_t L, size_t U>
inline HALMD_GPU_ENABLED
typename std::enable_if<(U <= L), void>::type
element_prod(fixed_vector<T, N>& v, fixed_vector<T, N> const& w)
{
    v[L] *= w[L];
}

template <typename T, size_t N, size_t L, size_t U>
inline HALMD_GPU_ENABLED
typename std::enable_if<(U > L), void>::type
element_prod(fixed_vector<T, N>& v, fixed_vector<T, N> const& w)
{
    element_prod<T, N, L, (L + U) / 2>(v, w);
    element_prod<T, N, (L + U) / 2 + 1, U>(v, w);
}

template <typename T, size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<T, N> element_prod(fixed_vector<T, N> v, fixed_vector<T, N> const& w)
{
    element_prod<T, N, 0, N - 1>(v, w);
    return v;
}

/**
 * Elementwise vector division
 */
template <typename T, size_t N, size_t L, size_t U>
inline HALMD_GPU_ENABLED
typename std::enable_if<(U <= L), void>::type
element_div(fixed_vector<T, N>& v, fixed_vector<T, N> const& w)
{
    v[L] /= w[L];
}

template <typename T, size_t N, size_t L, size_t U>
inline HALMD_GPU_ENABLED
typename std::enable_if<(U > L), void>::type
element_div(fixed_vector<T, N>& v, fixed_vector<T, N> const& w)
{
    element_div<T, N, L, (L + U) / 2>(v, w);
    element_div<T, N, (L + U) / 2 + 1, U>(v, w);
}

template <typename T, size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<T, N> element_div(fixed_vector<T, N> v, fixed_vector<T, N> const& w)
{
    element_div<T, N, 0, N - 1>(v, w);
    return v;
}

/**
 * Elementwise vector modulus
 */
template <typename T, size_t N, size_t L, size_t U>
inline HALMD_GPU_ENABLED
typename std::enable_if<(U <= L), void>::type
element_mod(fixed_vector<T, N>& v, fixed_vector<T, N> const& w)
{
    v[L] %= w[L];
}

template <typename T, size_t N, size_t L, size_t U>
inline HALMD_GPU_ENABLED
typename std::enable_if<(U > L), void>::type
element_mod(fixed_vector<T, N>& v, fixed_vector<T, N> const& w)
{
    element_mod<T, N, L, (L + U) / 2>(v, w);
    element_mod<T, N, (L + U) / 2 + 1, U>(v, w);
}

template <typename T, size_t N>
inline HALMD_GPU_ENABLED
typename std::enable_if<std::is_integral<T>::value, fixed_vector<T, N> >::type
element_mod(fixed_vector<T, N> v, fixed_vector<T, N> const& w)
{
    element_mod<T, N, 0, N - 1>(v, w);
    return v;
}

/**
 * Elementwise minimum
 */
template <typename T, size_t N, size_t L, size_t U>
inline HALMD_GPU_ENABLED
typename std::enable_if<(U <= L), void>::type
element_min(fixed_vector<T, N>& v, fixed_vector<T, N> const& w)
{
    HALMD_GPU_USING(::min, std::min);
    v[L] = min(v[L], w[L]);
}

template <typename T, size_t N, size_t L, size_t U>
inline HALMD_GPU_ENABLED
typename std::enable_if<(U > L), void>::type
element_min(fixed_vector<T, N>& v, fixed_vector<T, N> const& w)
{
    element_min<T, N, L, (L + U) / 2>(v, w);
    element_min<T, N, (L + U) / 2 + 1, U>(v, w);
}

template <typename T, size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<T, N> element_min(fixed_vector<T, N> v, fixed_vector<T, N> const& w)
{
    element_min<T, N, 0, N - 1>(v, w);
    return v;
}

/**
 * Elementwise maximum
 */
template <typename T, size_t N, size_t L, size_t U>
inline HALMD_GPU_ENABLED
typename std::enable_if<(U <= L), void>::type
element_max(fixed_vector<T, N>& v, fixed_vector<T, N> const& w)
{
    HALMD_GPU_USING(::max, std::max);
    v[L] = max(v[L], w[L]);
}

template <typename T, size_t N, size_t L, size_t U>
inline HALMD_GPU_ENABLED
typename std::enable_if<(U > L), void>::type
element_max(fixed_vector<T, N>& v, fixed_vector<T, N> const& w)
{
    element_max<T, N, L, (L + U) / 2>(v, w);
    element_max<T, N, (L + U) / 2 + 1, U>(v, w);
}

template <typename T, size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<T, N> element_max(fixed_vector<T, N> v, fixed_vector<T, N> const& w)
{
    element_max<T, N, 0, N - 1>(v, w);
    return v;
}

/**
 * Elementwise round to nearest integer not greater than argument
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED
typename std::enable_if<std::is_same<T, double>::value, fixed_vector<T, N> >::type
floor(fixed_vector<T, N> v)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] = ::floor(v[i]);
    }
    return v;
}

/**
 * Elementwise round to nearest integer not greater than argument
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED
typename std::enable_if<std::is_same<T, float>::value, fixed_vector<T, N> >::type
floor(fixed_vector<T, N> v)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] = ::floorf(v[i]);
    }
    return v;
}

/**
 * Elementwise round to nearest integer not less argument
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED
typename std::enable_if<std::is_same<T, double>::value, fixed_vector<T, N> >::type
ceil(fixed_vector<T, N> v)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] = ::ceil(v[i]);
    }
    return v;
}

/**
 * Elementwise round to nearest integer not less argument
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED
typename std::enable_if<std::is_same<T, float>::value, fixed_vector<T, N> >::type
ceil(fixed_vector<T, N> v)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] = ::ceilf(v[i]);
    }
    return v;
}

/**
 * Elementwise square root function
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED
typename std::enable_if<std::is_floating_point<T>::value, fixed_vector<T, N> >::type
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
typename std::enable_if<std::is_floating_point<T>::value, fixed_vector<T, N> >::type
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
typename std::enable_if<std::is_floating_point<T>::value, fixed_vector<T, N> >::type
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
typename std::enable_if<std::is_floating_point<T>::value, fixed_vector<T, N> >::type
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
typename std::enable_if<std::is_floating_point<T>::value, fixed_vector<T, N> >::type
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
typename std::enable_if<std::is_same<T, double>::value, fixed_vector<T, N> >::type
rint(fixed_vector<T, N> v)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] = ::rint(v[i]);
    }
    return v;
}

template <typename T, size_t N>
inline HALMD_GPU_ENABLED
typename std::enable_if<std::is_same<T, float>::value, fixed_vector<T, N> >::type
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
typename std::enable_if<std::is_same<T, double>::value, fixed_vector<T, N> >::type
round(fixed_vector<T, N> v)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] = ::round(v[i]);
    }
    return v;
}

template <typename T, size_t N>
inline HALMD_GPU_ENABLED
typename std::enable_if<std::is_same<T, float>::value, fixed_vector<T, N> >::type
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

} // namespace detail
} // namespace blas
} // namespace numeric
} // namespace halmd

#endif /* ! HALMD_NUMERIC_BLAS_DETAIL_OPERATORS_HPP */
