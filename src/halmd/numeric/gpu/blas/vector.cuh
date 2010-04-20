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

#ifndef HALMD_NUMERIC_GPU_BLAS_VECTOR_CUH
#define HALMD_NUMERIC_GPU_BLAS_VECTOR_CUH

#include <halmd/numeric/gpu/blas/detail/vector2d.cuh>
#include <halmd/numeric/gpu/blas/detail/vector3d.cuh>
#include <halmd/numeric/gpu/blas/detail/vector4d.cuh>

namespace halmd { namespace numeric { namespace gpu { namespace blas
{

using detail::vector;

}}}} // namespace halmd::numeric::gpu::blas

#endif /* ! HALMD_NUMERIC_GPU_BLAS_VECTOR_CUH */
