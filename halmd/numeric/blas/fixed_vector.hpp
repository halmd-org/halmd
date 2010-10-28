/*
 * Copyright © 2008-2010  Peter Colberg and Felix Höfling
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

#include <halmd/numeric/blas/fixed_vector/size_N.hpp>
#include <halmd/numeric/blas/fixed_vector/size_2.hpp>
#include <halmd/numeric/blas/fixed_vector/size_3.hpp>
#include <halmd/numeric/blas/fixed_vector/size_4.hpp>
#include <halmd/numeric/blas/fixed_vector/size_6.hpp>
#include <halmd/numeric/blas/fixed_vector/operators.hpp>
#include <halmd/numeric/blas/fixed_vector/rounding.hpp>

namespace halmd
{

// import into top-level namespace
using detail::numeric::blas::fixed_vector;

} // namespace halmd

#endif /* ! HALMD_NUMERIC_BLAS_FIXED_VECTOR_HPP */
