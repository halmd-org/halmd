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

#ifndef HALMD_NUMERIC_BLAS_DETAIL_STORAGE_HPP
#define HALMD_NUMERIC_BLAS_DETAIL_STORAGE_HPP

#include <halmd/config.hpp>

#ifndef __CUDACC__
# include <boost/array.hpp>
#endif

namespace halmd {
namespace numeric {
namespace blas {
namespace detail {

#ifndef __CUDACC__

template <typename T, size_t N>
struct fixed_array
  : boost::array<T, N> {};

#else /* __CUDACC__ */

//
// The purpose of a fixed_array is to serve as the underlying
// array type to a fixed-length algebraic vector. It defines
// operator[] to allow convenient access of its components.
//
template <typename T, size_t N>
struct fixed_array
{
public:
    typedef T value_type;
    typedef value_type& reference;
    typedef value_type const& const_reference;
    typedef size_t size_type;
    enum { static_size = N };

    HALMD_GPU_ENABLED reference operator[](size_type i)
    {
        return storage_[i];
    }

    HALMD_GPU_ENABLED const_reference operator[](size_type i) const
    {
        return storage_[i];
    }

private:
    T storage_[N];
};

#endif /* __CUDACC__ */

} // namespace detail
} // namespace blas
} // namespace numeric
} // namespace halmd

#endif /* ! HALMD_NUMERIC_BLAS_DETAIL_STORAGE_HPP */
