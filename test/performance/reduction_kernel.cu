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

#include <halmd/algorithm/gpu/reduction.cuh>
#include <halmd/algorithm/gpu/transform.cuh>
#include <halmd/numeric/mp/dsfloat.hpp>
#include <halmd/utility/gpu/thread.cuh>
#include <test/performance/reduction_kernel.hpp>

#include <cuda_wrapper/cuda_wrapper.hpp>

using namespace halmd::algorithm::gpu;

template <typename T>
__global__ void reduce_sum(T const* input, T* output)
{
    for (int i = 0; i < NREDUCES; ++i) {
        T val = input[TID + i * NTHREADS];
        reduce<sum_>(val);

        if (TID == 0) {
            output[i] = val;
        }
    }
}

cuda::function<void (float const*, float*)> reduce_float_kernel(reduce_sum<float>);
cuda::function<void (int const*, int*)> reduce_int_kernel(reduce_sum<int>);
cuda::function<void (halmd::dsfloat const*, halmd::dsfloat*)> reduce_dsfloat_kernel(reduce_sum<halmd::dsfloat>);
cuda::function<void (halmd::fixed_vector<float, 3> const*, halmd::fixed_vector<float, 3>*)>
    reduce_fixed_vector_float_kernel(reduce_sum<halmd::fixed_vector<float, 3>>);
cuda::function<void (halmd::fixed_vector<halmd::dsfloat, 3> const*, halmd::fixed_vector<halmd::dsfloat, 3>*)>
    reduce_fixed_vector_dsfloat_kernel(reduce_sum<halmd::fixed_vector<halmd::dsfloat, 3>>);
