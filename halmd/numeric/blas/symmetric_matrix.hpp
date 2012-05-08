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

#ifndef HALMD_NUMERIC_BLAS_SYMMETRIC_MATRIX_HPP
#define HALMD_NUMERIC_BLAS_SYMMETRIC_MATRIX_HPP

#include <halmd/config.hpp>

namespace halmd {
namespace symmetric_matrix {

/**
 * Compute storage index for lower symmetric matrix
 */
inline HALMD_GPU_ENABLED unsigned int lower_index(unsigned int row, unsigned int col)
{
    unsigned int row_lower = col > row ? col : row;
    unsigned int col_lower = col < row ? col : row;
    return col_lower + ((row_lower * (row_lower + 1)) / 2);
}

} // namespace symmetric_matrix
} // namespace halmd

#endif /* ! HALMD_NUMERIC_BLAS_SYMMETRIC_MATRIX_HPP */
