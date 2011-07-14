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

#ifndef HALMD_NUMERIC_BLAS_FIXED_VECTOR_SIZE_N_HPP
#define HALMD_NUMERIC_BLAS_FIXED_VECTOR_SIZE_N_HPP

#include <boost/type_traits/is_convertible.hpp>
#include <boost/utility/enable_if.hpp>

#include <halmd/config.hpp>
#include <halmd/numeric/blas/fixed_array.hpp>

namespace halmd {
namespace detail {
namespace numeric {
namespace blas {

/**
 * N-dimensional vector of arbitrary value type
 */
template <typename T, size_t N>
struct fixed_vector
  : fixed_array<T, N>
{
    typedef fixed_array<T, N> _Base;
    typedef typename _Base::value_type value_type;
    enum { static_size = _Base::static_size };

    HALMD_GPU_ENABLED fixed_vector() {}

    /**
     * Initialization by scalar
     */
    template <typename U>
    HALMD_GPU_ENABLED fixed_vector(U const& s,
      typename boost::enable_if<boost::is_convertible<U, T> >::type* dummy = 0)
    {
        for (size_t i = 0; i < static_size; ++i) {
            (*this)[i] = s;
        }
    }

    /**
     * Explicit conversion from vector of convertible element type
     */
    template <typename U>
    HALMD_GPU_ENABLED explicit fixed_vector(fixed_vector<U, N> const& v,
      typename boost::enable_if<boost::is_convertible<U, T> >::type* dummy = 0)
    {
        for (size_t i = 0; i < static_size; ++i) {
            (*this)[i] = static_cast<value_type>(v[i]);
        }
    }
};

} // namespace detail
} // namespace numeric
} // namespace blas
// import into top-level namespace
using detail::numeric::blas::fixed_vector;

} // namespace halmd

#endif /* ! HALMD_NUMERIC_BLAS_FIXED_VECTOR_SIZE_N_HPP */
