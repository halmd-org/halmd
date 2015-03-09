/*
 * Copyright © 2015 Nicolas Höft
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

#ifndef HALMD_ALGORITHM_GPU_COPY_IF_KERNEL_CUH
#define HALMD_ALGORITHM_GPU_COPY_IF_KERNEL_CUH

#include <cub/device/device_select.cuh>
#include <cuda_wrapper/error.hpp>

#include <halmd/algorithm/gpu/copy_if_kernel.hpp>
#include <halmd/utility/gpu/caching_array.cuh>

namespace halmd {
namespace algorithm {
namespace gpu {
namespace copy_if_kernel {

template <typename InputIterator, typename OutputIterator, typename Predicate>
unsigned int copy_if(InputIterator g_input, unsigned int size, Predicate predicate, OutputIterator g_output)
{
    caching_array<int> g_num_selected_out(1);

    // Determine temporary device storage requirements
    size_t temp_storage_bytes = 0;
    CUDA_CALL(cub::DeviceSelect::If(
        0
      , temp_storage_bytes
      , g_input
      , g_output
      , g_num_selected_out.begin()
      , size
      , predicate
    ));

    caching_array<char> g_temp_storage(temp_storage_bytes);

    // copy the elements according to the predicate
    CUDA_CALL(cub::DeviceSelect::If(
        g_temp_storage.begin()
      , temp_storage_bytes
      , g_input
      , g_output
      , g_num_selected_out.begin()
      , size
      , predicate
    ));

    // the output length resides in gpu, copy it to host memory
    int h_output_len;
    CUDA_CALL(cudaMemcpy(
        &h_output_len
      , g_num_selected_out.begin()
      , sizeof(int)
      , cudaMemcpyDeviceToHost
    ));
    return h_output_len;
}

} // namespace copy_if_kernel

/**
 * CUDA C++ wrapper
 */
template <typename InputIterator, typename OutputIterator, typename Predicate>
copy_if_wrapper<InputIterator, OutputIterator, Predicate> const copy_if_wrapper<InputIterator, OutputIterator, Predicate>::kernel = {
   copy_if_kernel::copy_if<InputIterator, OutputIterator, Predicate>
};

} // namespace gpu
} // namespace algorithm
} // namespace halmd

#endif /* ! HALMD_ALGORITHM_GPU_COPY_IF_KERNEL_CUH */
