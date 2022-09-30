/*
 * Copyright © 2021 Jaslo Ziska
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

#ifndef HALMD_ALGORITHM_GPU_BITS_SHFL_CUH
#define HALMD_ALGORITHM_GPU_BITS_SHFL_CUH

#include <type_traits>

#include <halmd/numeric/blas/fixed_vector.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>

//
// Overloaded shfl functions for HALMD specific types
//

namespace halmd {
namespace algorithm {
namespace gpu {
namespace bits {

/**
 * Shuffle down (sync)
 *
 * See CUDA documentation for details:
 * <https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#warp-shuffle-functions>
 */
template <typename T>
__inline__ __device__ typename std::enable_if<std::is_integral<T>::value || std::is_floating_point<T>::value, T>::type
shfl_down(unsigned int mask, T var, unsigned int delta, int width = warpSize)
{
    return __shfl_down_sync(mask, var, delta, width);
}

/**
 * Shuffle for types which are not supported by __shfl_down_sync
 *
 * Inspired by the "ShuffleDown" function fron the cub library, file "cub/util_ptx.cuh", see:
 * <https://github.com/NVIDIA/cub/blob/main/cub/util_ptx.cuh#L550>
 */
template <typename T>
__inline__ __device__ typename std::enable_if<!(std::is_integral<T>::value || std::is_floating_point<T>::value), T>::type
shfl_down(unsigned int mask, T val, unsigned int delta, int width = warpSize)
{
    typedef typename std::conditional<
        (sizeof(T) % sizeof(unsigned long long) == 0 && alignof(T) % alignof(unsigned long long) == 0)
      , unsigned long long
      , typename std::conditional<
            (sizeof(T) % sizeof(unsigned int) == 0 && alignof(T) % alignof(unsigned int) == 0)
          , unsigned int
          , typename std::conditional<
                (sizeof(T) % sizeof(unsigned short) == 0 && alignof(T) % alignof(unsigned short) == 0)
              , unsigned short
              , unsigned char
            >::type
        >::type
    >::type shfl_type;

    shfl_type* parts = reinterpret_cast<shfl_type*>(&val);

    for (int i = 0; i < (sizeof(T) + sizeof(shfl_type) - 1) / sizeof(shfl_type); ++i) {
        parts[i] = __shfl_down_sync(mask, parts[i], delta, width);
    }

    return val;
}

/**
 * Overloaded shuffle for halmd::dsfloat
 */
__inline__ __device__ halmd::dsfloat
shfl_down(unsigned mask, halmd::dsfloat var, unsigned int delta, int width = warpSize)
{
    var.hi = __shfl_down_sync(mask, var.hi, delta, width);
    var.lo = __shfl_down_sync(mask, var.lo, delta, width);

    return var;
}

/**
 * Overloaded shuffle for halmd::fixed_vector
 */
template <typename T, size_t N>
__inline__ __device__ halmd::fixed_vector<T, N>
shfl_down(unsigned mask, halmd::fixed_vector<T, N> var, unsigned int delta, int width = warpSize)
{
    for (int i = 0; i < N; ++i) {
        var[i] = shfl_down(mask, var[i], delta, width);
    }

    return var;
}

} // namespace bits
} // namespace gpu
} // namespace algorithm
} // namespace halmd

#endif // ! HALMD_ALGORITHM_GPU_BITS_SHFL_CUH
