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

#ifndef HALMD_NUMERIC_HOST_BLAS_DETAIL_BLAS1_HPP
#define HALMD_NUMERIC_HOST_BLAS_DETAIL_BLAS1_HPP

namespace halmd { namespace numeric { namespace host { namespace blas
{

namespace detail
{

/** provide functionality of BLAS Level 1 */

// see e.g. boost/numeric/ublas/blas.hpp

/**
 * absolute sum or 1-norm, sum(|x_i|)
 */
template <typename vector_type>
inline typename vector_type::value_type norm_1(vector_type const& v)
{
    typename vector_type::value_type s = 0;
    for (typename vector_type::size_type i=0; i < v.size(); i++) {
        s += std::abs(v[i]);
    }
    return s;
}

/**
 * Euclidean norm or 2-norm, sqrt(sum(x_i^2))
 */
template <typename vector_type>
inline typename vector_type::value_type norm_2(vector_type const& v)
{
    return std::sqrt(inner_prod(v, v));
}

/**
 * maximum norm or inf-norm, max(|x_i|)
 */
template <typename vector_type>
inline typename vector_type::value_type norm_inf(vector_type const& v)
{
    typename vector_type::value_type s = 0;
    for (typename vector_type::size_type i=0; i < v.size(); i++) {
        s = std::max(s, std::abs(v[i]));
    }
    return s;
}

/**
 * index of maximal element, first i where |x_i| = max(|x_i|)
 */
template <typename vector_type>
inline typename vector_type::size_type index_norm_inf(vector_type const& v)
{
    typedef typename vector_type::value_type value_type;
    typedef typename vector_type::size_type size_type;
    value_type s = 0;
    size_type k = 0;
    for (size_type i=0; i < v.size(); i++) {
        value_type t = std::abs(v[i]);
        if (s > t) {
            t = s;
            k = i;
        }
    }
    return k;
}

/** TODO: implement plane rotations */

} // namespace detail

}}}} // namespace halmd::numeric::host::blas

#endif /* ! HALMD_NUMERIC_HOST_BLAS_DETAIL_BLAS1_HPP */
