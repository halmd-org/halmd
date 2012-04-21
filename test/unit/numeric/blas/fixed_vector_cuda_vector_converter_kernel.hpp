/*
 * Copyright © 2012  Peter Colberg
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

#ifndef TEST_UNIT_NUMERIC_BLAS_FIXED_VECTOR_CUDA_VECTOR_CONVERTER_KERNEL_HPP
#define TEST_UNIT_NUMERIC_BLAS_FIXED_VECTOR_CUDA_VECTOR_CONVERTER_KERNEL_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>

template <typename U, typename V>
struct float_kernel
{
    cuda::function<void (float4 const*, float4*, U*, V*) > converter_one;
    cuda::function<void (float4 const*, float4*, U*, V*) > converter_two;
    static float_kernel const kernel;
};

template <typename U, typename V>
struct double_kernel
{
    cuda::function<void (float4 const*, float4*, U*, V*) > converter_two;
    static double_kernel const kernel;
};

#endif /* ! TEST_UNIT_NUMERIC_BLAS_FIXED_VECTOR_CUDA_VECTOR_CONVERTER_KERNEL_HPP */
