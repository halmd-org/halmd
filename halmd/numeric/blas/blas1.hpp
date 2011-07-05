/*
 * Copyright © 2010  Felix Höfling
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

#ifndef HALMD_NUMERIC_BLAS_BLAS1_HPP
#define HALMD_NUMERIC_BLAS_BLAS1_HPP

#include <halmd/config.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>

namespace halmd {
namespace detail {
namespace numeric {
namespace blas {

/** provide functionality of BLAS Level 1 */

// see e.g. boost/numeric/ublas/blas.hpp

/**
 * absolute sum or 1-norm, sum(|x_i|)
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED T norm_1(fixed_vector<T, N> const& v)
{
    HALMD_GPU_USING(::abs, std::abs);
    T s = 0;
    for (size_t i = 0; i < v.size(); ++i) {
        s += abs(v[i]);
    }
    return s;
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
template <typename T, size_t N>
inline HALMD_GPU_ENABLED T norm_inf(fixed_vector<T, N> const& v)
{
    HALMD_GPU_USING(::abs, std::abs);
    HALMD_GPU_USING(::max, std::max);
    T s = 0;
    for (size_t i = 0; i < v.size(); ++i) {
        s = max(s, abs(v[i]));
    }
    return s;
}

/**
 * index of maximal element, first i where |x_i| = max(|x_i|)
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED size_t index_norm_inf(fixed_vector<T, N> const& v)
{
    HALMD_GPU_USING(::abs, std::abs);
    T s = 0;
    size_t k = 0;
    for (size_t i = 0; i < v.size(); ++i) {
        T t = abs(v[i]);
        if (s > t) {
            t = s;
            k = i;
        }
    }
    return k;
}

/** TODO: implement plane rotations */

} // namespace detail
} // namespace numeric
} // namespace blas
} // namespace halmd

#endif /* ! HALMD_NUMERIC_BLAS_BLAS1_HPP */
