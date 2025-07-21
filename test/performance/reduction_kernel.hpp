/*
 * Copyright © 2021 Jaslo Ziska
 * Copyright © 2023 Felix Höfling
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

#ifndef HALMD_TEST_PERFORMANCE_REDUCTION_KERNEL_HPP
#define HALMD_TEST_PERFORMANCE_REDUCTION_KERNEL_HPP

#include <halmd/numeric/blas/fixed_vector.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>

#include <cuda_wrapper/cuda_wrapper.hpp>

static const size_t NTHREADS = 224; // NTHREADS / 32 is not a power of two
static const size_t NREDUCES = 5000;

template <typename T>
struct reduce_kernel
{
    cuda::function<void (T const*, T*)> sum;
    static reduce_kernel kernel;
};

#endif // ! HALMD_TEST_PERFORMANCE_REDUCTION_KERNEL_HPP
