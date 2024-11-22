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

#ifndef HALMD_NUMERIC_BLAS_DETAIL_VECTOR_HPP
#define HALMD_NUMERIC_BLAS_DETAIL_VECTOR_HPP

#include <halmd/config.hpp>

#include <algorithm>
#include <cassert>
#include <initializer_list>
#include <type_traits>

#include <halmd/numeric/blas/detail/array.hpp>

namespace halmd {
namespace numeric {
namespace blas {
namespace detail {

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
     * Assign values from initializer list.
     */
    template <typename U>
    HALMD_GPU_ENABLED fixed_vector(std::initializer_list<U> const& v,
      typename std::enable_if<std::is_convertible<U, T>::value>::type* dummy = 0)
    {
        assert( v.size() == _Base::size() );
        std::copy(v.begin(), v.end(), _Base::begin());
    }

    /**
     * Initialization by scalar
     */
    HALMD_GPU_ENABLED fixed_vector(T const& s)
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
      typename std::enable_if<std::is_convertible<U, T>::value>::type* dummy = 0)
    {
        for (size_t i = 0; i < static_size; ++i) {
            (*this)[i] = static_cast<value_type>(v[i]);
        }
    }
};

} // namespace detail
} // namespace blas
} // namespace numeric
} // namespace halmd

#endif /* ! HALMD_NUMERIC_BLAS_DETAIL_VECTOR_HPP */
