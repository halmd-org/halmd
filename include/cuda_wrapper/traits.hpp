/*
 * Copyright (C) 2010  Peter Colberg
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

#ifndef CUDA_TRAITS_HPP
#define CUDA_TRAITS_HPP

#include <cstddef>
#include <cuda.h>

namespace cuda
{

#if (CUDA_VERSION >= 3020)
typedef std::size_t size_type;
#else
typedef unsigned int size_type;
#endif

} // namespace cuda

#endif /* ! CUDA_TRAITS_HPP */
