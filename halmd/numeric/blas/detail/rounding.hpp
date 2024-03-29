/*
 * Copyright © 2008-2010  Peter Colberg
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

#ifndef HALMD_NUMERIC_BLAS_DETAIL_ROUNDING_HPP
#define HALMD_NUMERIC_BLAS_DETAIL_ROUNDING_HPP

#include <halmd/config.hpp>

#include <type_traits>
#ifndef __CUDACC__
# include <cmath>
#endif
#ifdef HALMD_WITH_GPU
# include <cuda_runtime.h> // CUDA vector types for host compiler
#endif

#include <halmd/numeric/blas/detail/vector.hpp>

namespace halmd {
namespace numeric {
namespace blas {
namespace detail {

#ifdef __CUDACC__

/**
 * Fast, accurate floating-point division by s < 2^126
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED
typename std::enable_if<std::is_same<T, double>::value, fixed_vector<T, N>>::type
fdivide(fixed_vector<T, N> v, fixed_vector<T, N> const& w)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] = ::fdivide(v[i], w[i]);
    }
    return v;
}

template <typename T, size_t N>
inline HALMD_GPU_ENABLED
typename std::enable_if<std::is_same<T, float>::value, fixed_vector<T, N>>::type
fdivide(fixed_vector<T, N> v, fixed_vector<T, N> const& w)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] = ::fdividef(v[i], w[i]);
    }
    return v;
}

/**
 * Limit floating-point elements to unit interval [0, 1]
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED
typename std::enable_if<std::is_same<T, float>::value, fixed_vector<T, N>>::type
saturate(fixed_vector<T, N> v)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] = ::saturate(v[i]);
    }
    return v;
}

/**
 * Convert floating-point elements to integers, rounding to negative infinity
 */
#ifdef HALMD_GPU_DOUBLE_PRECISION
template <size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<int, N> __double2int_rd(fixed_vector<double, N> const& v)
{
    fixed_vector<int, N> w;
    for (size_t i = 0; i < N; ++i) {
        w[i] = ::__double2int_rd(v[i]);
    }
    return w;
}
#endif

template <size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<int, N> __float2int_rd(fixed_vector<float, N> const& v)
{
    fixed_vector<int, N> w;
    for (size_t i = 0; i < N; ++i) {
        w[i] = ::__float2int_rd(v[i]);
    }
    return w;
}

/**
 * Convert floating-point elements to integers, rounding to nearest even integer
 */
#ifdef HALMD_GPU_DOUBLE_PRECISION
template <size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<int, N> __double2int_rn(fixed_vector<double, N> const& v)
{
    fixed_vector<int, N> w;
    for (size_t i = 0; i < N; ++i) {
        w[i] = ::__double2int_rn(v[i]);
    }
    return w;
}
#endif

template <size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<int, N> __float2int_rn(fixed_vector<float, N> const& v)
{
    fixed_vector<int, N> w;
    for (size_t i = 0; i < N; ++i) {
        w[i] = ::__float2int_rn(v[i]);
    }
    return w;
}

/**
 * Convert floating-point elements to integers, rounding to positive infinity
 */
#ifdef HALMD_GPU_DOUBLE_PRECISION
template <size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<int, N> __double2int_ru(fixed_vector<double, N> const& v)
{
    fixed_vector<int, N> w;
    for (size_t i = 0; i < N; ++i) {
        w[i] = ::__double2int_ru(v[i]);
    }
    return w;
}
#endif

template <size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<int, N> __float2int_ru(fixed_vector<float, N> const& v)
{
    fixed_vector<int, N> w;
    for (size_t i = 0; i < N; ++i) {
        w[i] = ::__float2int_ru(v[i]);
    }
    return w;
}

/**
 * Convert floating-point elements to integers, rounding towards zero
 */
#ifdef HALMD_GPU_DOUBLE_PRECISION
template <size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<int, N> __double2int_rz(fixed_vector<double, N> const& v)
{
    fixed_vector<int, N> w;
    for (size_t i = 0; i < N; ++i) {
        w[i] = ::__double2int_rz(v[i]);
    }
    return w;
}
#endif

template <size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<int, N> __float2int_rz(fixed_vector<float, N> const& v)
{
    fixed_vector<int, N> w;
    for (size_t i = 0; i < N; ++i) {
        w[i] = ::__float2int_rz(v[i]);
    }
    return w;
}

/**
 * Convert floating-point elements to unsigned integers, rounding to negative infinity
 */
#ifdef HALMD_GPU_DOUBLE_PRECISION
template <size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<unsigned int, N> __double2uint_rd(fixed_vector<double, N> const& v)
{
    fixed_vector<unsigned int, N> w;
    for (size_t i = 0; i < N; ++i) {
        w[i] = ::__double2uint_rd(v[i]);
    }
    return w;
}
#endif

template <size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<unsigned int, N> __float2uint_rd(fixed_vector<float, N> const& v)
{
    fixed_vector<unsigned int, N> w;
    for (size_t i = 0; i < N; ++i) {
        w[i] = ::__float2uint_rd(v[i]);
    }
    return w;
}

/**
 * Convert floating-point elements to unsigned integers, rounding to nearest even integer
 */
#ifdef HALMD_GPU_DOUBLE_PRECISION
template <size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<unsigned int, N> __double2uint_rn(fixed_vector<double, N> const& v)
{
    fixed_vector<unsigned int, N> w;
    for (size_t i = 0; i < N; ++i) {
        w[i] = ::__double2uint_rn(v[i]);
    }
    return w;
}
#endif

template <size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<unsigned int, N> __float2uint_rn(fixed_vector<float, N> const& v)
{
    fixed_vector<unsigned int, N> w;
    for (size_t i = 0; i < N; ++i) {
        w[i] = ::__float2uint_rn(v[i]);
    }
    return w;
}

/**
 * Convert floating-point elements to unsigned integers, rounding to positive infinity
 */
#ifdef HALMD_GPU_DOUBLE_PRECISION
template <size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<unsigned int, N> __double2uint_ru(fixed_vector<double, N> const& v)
{
    fixed_vector<unsigned int, N> w;
    for (size_t i = 0; i < N; ++i) {
        w[i] = ::__double2uint_ru(v[i]);
    }
    return w;
}
#endif

template <size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<unsigned int, N> __float2uint_ru(fixed_vector<float, N> const& v)
{
    fixed_vector<unsigned int, N> w;
    for (size_t i = 0; i < N; ++i) {
        w[i] = ::__float2uint_ru(v[i]);
    }
    return w;
}

/**
 * Convert floating-point elements to unsigned integers, rounding towards zero
 */
#ifdef HALMD_GPU_DOUBLE_PRECISION
template <size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<unsigned int, N> __double2uint_rz(fixed_vector<double, N> const& v)
{
    fixed_vector<unsigned int, N> w;
    for (size_t i = 0; i < N; ++i) {
        w[i] = ::__double2uint_rz(v[i]);
    }
    return w;
}
#endif

template <size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<unsigned int, N> __float2uint_rz(fixed_vector<float, N> const& v)
{
    fixed_vector<unsigned int, N> w;
    for (size_t i = 0; i < N; ++i) {
        w[i] = ::__float2uint_rz(v[i]);
    }
    return w;
}

#endif /* ! __CUDACC__ */

} // namespace detail
} // namespace blas
} // namespace numeric
} // namespace halmd

#endif /* ! HALMD_NUMERIC_BLAS_DETAIL_ROUNDING_HPP */
