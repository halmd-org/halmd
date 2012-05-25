/*
 * Copyright © 2010  Felix Höfling and Peter Colberg
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

#ifndef HALMD_IO_UTILITY_HDF5_HPP
#define HALMD_IO_UTILITY_HDF5_HPP

#include <boost/numeric/ublas/storage.hpp>

#include <h5xx/h5xx.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>
#include <halmd/utility/raw_array.hpp>

namespace h5xx {

template <typename T, size_t N>
struct is_array<halmd::numeric::blas::detail::fixed_vector<T, N> >
{
    enum { value = true };
    typedef is_array type;
};

template <typename T>
struct is_vector<halmd::raw_array<T> >
{
    enum { value = true };
    typedef is_vector type;
};

template <typename T, typename Alloc>
struct is_vector<boost::numeric::ublas::unbounded_array<T, Alloc> >
{
    enum { value = true };
    typedef is_vector type;
};

} // namespace h5xx

#endif /* ! HALMD_IO_UTILITY_HDF5_HPP */
