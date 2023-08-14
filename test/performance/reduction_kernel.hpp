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

#ifndef HALMD_TEST_PERFORMANCE_REDUCTION_KERNEL_HPP
#define HALMD_TEST_PERFORMANCE_REDUCTION_KERNEL_HPP

#include <halmd/numeric/blas/fixed_vector.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>

#include <cuda_wrapper/cuda_wrapper.hpp>

static const size_t NTHREADS = 224; // NTHREADS / 32 is not a power of two
static const size_t NREDUCES = 10000;

extern cuda::function<void (float const*, float*)> reduce_float_kernel;
extern cuda::function<void (int const*, int*)> reduce_int_kernel;
extern cuda::function<void (halmd::dsfloat const*, halmd::dsfloat*)> reduce_dsfloat_kernel;
extern cuda::function<void (halmd::fixed_vector<float, 3> const*, halmd::fixed_vector<float, 3>*)>
    reduce_fixed_vector_float_kernel;
extern cuda::function<void (halmd::fixed_vector<halmd::dsfloat, 3> const*, halmd::fixed_vector<halmd::dsfloat, 3>*)>
    reduce_fixed_vector_dsfloat_kernel;

#endif // ! HALMD_TEST_PERFORMANCE_REDUCTION_KERNEL_HPP
