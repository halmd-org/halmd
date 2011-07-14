/*
 * Copyright © 2008-2009  Peter Colberg
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

#ifndef HALMD_ALGORITHM_GPU_BITS_CUH
#define HALMD_ALGORITHM_GPU_BITS_CUH

//
// Bits and pieces used in GPU algorithms
//

namespace halmd {
namespace algorithm {
namespace gpu {

/**
 * swap arguments
 */
template <typename T>
__device__ __host__ void swap(T& a, T& b)
{
    T c = b;
    b = a;
    a = c;
}


} // namespace algorithm
} // namespace gpu
} // namespace halmd

#endif /* ! HALMD_ALGORITHM_GPU_BITS_CUH */
