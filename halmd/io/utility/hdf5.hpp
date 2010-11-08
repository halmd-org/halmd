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

#include <h5xx/h5xx.hpp>

namespace halmd { namespace detail { namespace numeric { namespace blas
{

// forward declaration
template <typename T, size_t N>
struct fixed_vector;

}}}} // namespace halmd::detail::numeric::blas

namespace h5xx
{

template <typename T, size_t N>
struct is_array<halmd::detail::numeric::blas::fixed_vector<T, N> >
  : boost::true_type {};

} // namespace h5xx

#endif /* ! HALMD_IO_UTILITY_HDF5_HPP */
