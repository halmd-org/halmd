/*
 * Copyright © 2011  Peter Colberg
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

#ifndef HALMD_TEST_UTILITY_GPU_VARIANT_HPP
#define HALMD_TEST_UTILITY_GPU_VARIANT_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>

/**
 * @file variant_kernel.hpp
 *
 * This file defines custom structs and C unions for testing on host
 * and GPU, and wraps the CUDA test kernels.
 */

//
// Data structures
//

/**
 * 4-component vector with default alignment
 */
template <typename T>
struct vector4
{
    T x, y, z, w;
};

/**
 * 4-component vector with 8-byte alignment
 *
 * Attribute packed needed if alignment of T is larger than 8 bytes,
 * e.g. for the test case with 16-byte aligned float4.
 */
template <typename T>
struct __attribute__((aligned(8), packed)) vector4_align8
{
    T x, y, z, w;
};

/**
 * 4-component vector with 16-byte alignment
 */
template <typename T>
struct __attribute__((aligned(16), packed)) vector4_align16
{
    T x, y, z, w;
};

//
// CUDA kernel wrappers
//

template <typename T>
struct prerequisites_wrapper
{
    cuda::function <void (size_t*)> alignment_of;
    cuda::function <void (size_t*)> sizeof_;
    static prerequisites_wrapper<T> kernel;
};

#endif /* !HALMD_TEST_UTILITY_GPU_VARIANT_HPP */
