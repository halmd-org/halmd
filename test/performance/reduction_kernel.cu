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

// initialise static class member
template <typename T>
reduce_kernel<T> reduce_kernel<T>::kernel = {
    reduce_sum<T>
};

// explicit template instantiations
using halmd::dsfloat;
using halmd::fixed_vector;

template class reduce_kernel<int>;
template class reduce_kernel<fixed_vector<int, 2>>;
#ifdef USE_GPU_SINGLE_PRECISION
template class reduce_kernel<float>;
template class reduce_kernel<fixed_vector<float, 3>>;
#endif
#ifdef USE_GPU_DOUBLE_SINGLE_PRECISION
template class reduce_kernel<dsfloat>;
template class reduce_kernel<fixed_vector<dsfloat, 2>>;
#endif
#ifdef USE_GPU_DOUBLE_PRECISION
template class reduce_kernel<double>;
template class reduce_kernel<fixed_vector<double, 3>>;
#endif

