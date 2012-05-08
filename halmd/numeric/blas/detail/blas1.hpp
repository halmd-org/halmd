/*
 * Copyright © 2010-2011  Felix Höfling and Peter Colberg
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

#ifndef HALMD_NUMERIC_BLAS_DETAIL_BLAS1_HPP
#define HALMD_NUMERIC_BLAS_DETAIL_BLAS1_HPP

#include <boost/mpl/int.hpp>
#include <boost/mpl/greater.hpp>
#include <boost/utility/enable_if.hpp>

#include <halmd/config.hpp>
#include <halmd/numeric/blas/detail/vector.hpp>

namespace halmd {
namespace numeric {
namespace blas {
namespace detail {

/** provide functionality of BLAS Level 1 */

// see e.g. boost/numeric/ublas/blas.hpp

/**
 * absolute sum or 1-norm, sum(|x_i|)
 */
template <typename T, size_t N, size_t L, size_t U>
inline HALMD_GPU_ENABLED
typename boost::disable_if<boost::mpl::greater<boost::mpl::int_<U>, boost::mpl::int_<L> >, T>::type
norm_1(fixed_vector<T, N> const& v)
{
    HALMD_GPU_USING(::abs, std::abs);
    return abs(v[L]);
}

template <typename T, size_t N, size_t L, size_t U>
inline HALMD_GPU_ENABLED
typename boost::enable_if<boost::mpl::greater<boost::mpl::int_<U>, boost::mpl::int_<L> >, T>::type
norm_1(fixed_vector<T, N> const& v)
{
    return norm_1<T, N, L, (L + U) / 2>(v) + norm_1<T, N, (L + U) / 2 + 1, U>(v);
}

template <typename T, size_t N>
inline HALMD_GPU_ENABLED T norm_1(fixed_vector<T, N> const& v)
{
    return norm_1<T, N, 0, N - 1>(v);
}

/**
 * Euclidean norm or 2-norm, sqrt(sum(x_i^2))
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED T norm_2(fixed_vector<T, N> const& v)
{
    HALMD_GPU_USING(::sqrt, std::sqrt);
    return sqrt(inner_prod(v, v));
}

/**
 * maximum norm or inf-norm, max(|x_i|)
 */
template <typename T, size_t N, size_t L, size_t U>
inline HALMD_GPU_ENABLED
typename boost::disable_if<boost::mpl::greater<boost::mpl::int_<U>, boost::mpl::int_<L> >, T>::type
norm_inf(fixed_vector<T, N> const& v)
{
    HALMD_GPU_USING(::abs, std::abs);
    return abs(v[L]);
}

template <typename T, size_t N, size_t L, size_t U>
inline HALMD_GPU_ENABLED
typename boost::enable_if<boost::mpl::greater<boost::mpl::int_<U>, boost::mpl::int_<L> >, T>::type
norm_inf(fixed_vector<T, N> const& v)
{
    HALMD_GPU_USING(::max, std::max);
    return max(norm_inf<T, N, L, (L + U) / 2>(v), norm_inf<T, N, (L + U) / 2 + 1, U>(v));
}

template <typename T, size_t N>
inline HALMD_GPU_ENABLED T norm_inf(fixed_vector<T, N> const& v)
{
    return norm_inf<T, N, 0, N - 1>(v);
}

/**
 * index of maximal element, first i where |x_i| = max(|x_i|)
 */
template <typename T, size_t N, size_t L, size_t U>
inline HALMD_GPU_ENABLED
typename boost::disable_if<boost::mpl::greater<boost::mpl::int_<U>, boost::mpl::int_<L> >, T>::type
index_norm_inf(fixed_vector<T, N> const& v, size_t& i)
{
    HALMD_GPU_USING(::abs, std::abs);
    i = L;
    return abs(v[L]);
}

template <typename T, size_t N, size_t L, size_t U>
inline HALMD_GPU_ENABLED
typename boost::enable_if<boost::mpl::greater<boost::mpl::int_<U>, boost::mpl::int_<L> >, T>::type
index_norm_inf(fixed_vector<T, N> const& v, size_t& i)
{
    size_t j, k;
    T s = index_norm_inf<T, N, L, (L + U) / 2>(v, j);
    T t = index_norm_inf<T, N, (L + U) / 2 + 1, U>(v, k);
    return t > s ? ((i = k), t) : ((i = j), s);
}

template <typename T, size_t N>
inline HALMD_GPU_ENABLED size_t index_norm_inf(fixed_vector<T, N> const& v)
{
    size_t i;
    index_norm_inf<T, N, 0, N - 1>(v, i);
    return i;
}

/** TODO: implement plane rotations */

} // namespace detail
} // namespace blas
} // namespace numeric
} // namespace halmd

#endif /* ! HALMD_NUMERIC_BLAS_DETAIL_BLAS1_HPP */
