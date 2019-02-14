/*
 * Copyright © 2012  Peter Colberg
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

#ifndef HALMD_UTILITY_GPU_SHARED_MEMORY_HPP
#define HALMD_UTILITY_GPU_SHARED_MEMORY_HPP

#include <boost/utility/enable_if.hpp>

#include <halmd/config.hpp> // HALMD_GPU_ARCH

namespace halmd {

/**
 * Maximum number of threads per block.
 */
struct max_threads_per_block
{
#if HALMD_GPU_ARCH >= 200
    enum { value = 1024 };
#else
    enum { value = 512 };
#endif
};

/**
 * Minimum number of threads per block.
 *
 * This is equivalent to the warp size.
 */
struct min_threads_per_block
{
    enum { value = 32 };
};

/**
 * Maximum shared memory in bytes
 */
struct max_shared_memory
{
#if HALMD_GPU_ARCH >= 200
    enum { value = 49152 };
#else
    enum { value = 16384 };
#endif
};

/**
 * Given an upper bound on the number of threads, calculates next-lower
 * power of 2. If the number of threads exceeds the maximum number of
 * threads per block, returns the maximum number of threads per block.
 */
template <unsigned int upper, unsigned int accumulator = 1, typename Enable = void>
struct power_of_two_max_threads;

template <unsigned int upper, unsigned int accumulator>
struct power_of_two_max_threads<upper, accumulator
  , typename boost::enable_if_c<
        (((upper > 1) && (accumulator == max_threads_per_block::value)) || (upper == 1))
    >::type>
{
    enum { value = accumulator };
};

template <unsigned int upper, unsigned int accumulator>
struct power_of_two_max_threads<upper, accumulator
  , typename boost::enable_if_c<
        ((upper > 1) && (accumulator < max_threads_per_block::value))
    >::type>
{
    enum { value = power_of_two_max_threads<(upper >> 1), (accumulator << 1)>::value };
};

/**
 * Given the element type of a shared memory array, calculates the maximum
 * number of threads per block supported by the GPU architecture, as the
 * next-lower power of 2.
 *
 * This metafunction calculates the maximum number of threads by dividing
 * the maximum shared memory by the element size, and then rounding the
 * result down to the next-lower power of 2.
 *
 * Note that with PTX version 1.3 and older, kernel arguments are passed
 * via shared memory. To account for the decrease in available shared
 * memory, this metafunction substracts 1 from the maximum number of
 * threads before rounding to the next-lower power of 2.
 * As an example, if the unrounded number of threads is exactly 512,
 * subtracting 1 will yield a next-lower power of two of 256 threads.
 *
 * For PTX version 2.0 and higher, kernel arguments are passed
 * via constant memory, thus this correction is not necessary.
 */
template <typename T>
struct shared_memory_max_threads
{
#if HALMD_GPU_ARCH >= 200
    enum { value = power_of_two_max_threads<max_shared_memory::value / sizeof(T)>::value };
#else
    enum { value = power_of_two_max_threads<(max_shared_memory::value / sizeof(T)) - 1>::value };
#endif
};

} // namespace halmd

#endif /* ! HALMD_UTILITY_GPU_SHARED_MEMORY_HPP */
