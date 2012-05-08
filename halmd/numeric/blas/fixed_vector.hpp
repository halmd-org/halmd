/*
 * Copyright © 2008-2012  Peter Colberg and Felix Höfling
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

#ifndef HALMD_NUMERIC_BLAS_FIXED_VECTOR_HPP
#define HALMD_NUMERIC_BLAS_FIXED_VECTOR_HPP

#include <halmd/config.hpp>

#include <halmd/numeric/blas/detail/vector.hpp>
#include <halmd/numeric/blas/detail/size_2.hpp>
#include <halmd/numeric/blas/detail/size_3.hpp>
#include <halmd/numeric/blas/detail/size_4.hpp>
#include <halmd/numeric/blas/detail/size_6.hpp>
#include <halmd/numeric/blas/detail/operators.hpp>
#include <halmd/numeric/blas/detail/rounding.hpp>
#ifdef CUDART_VERSION
# include <halmd/numeric/blas/detail/cuda_vector_converter.hpp>
#endif

namespace halmd {

#ifndef HALMD_NO_CXX11
// import into top-level namespace
template <typename T, size_t N>
using fixed_vector = halmd::numeric::blas::detail::fixed_vector<T, N>;
#else
using halmd::numeric::blas::detail::fixed_vector;
#endif

} // namespace halmd

#endif /* ! HALMD_NUMERIC_BLAS_FIXED_VECTOR_HPP */
