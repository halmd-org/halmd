/*
 * Copyright © 2011  Michael Kopp
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

#ifndef HALMD_ALGORITHM_GPU_APPLY_KERNEL_HPP
#define HALMD_ALGORITHM_GPU_APPLY_KERNEL_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>

#include <halmd/algorithm/gpu/transform.cuh>

namespace halmd {
namespace algorithm {
namespace gpu {

template <
    typename functor
  , typename input_type
  , typename coalesced_input_type       = input_type
  , typename output_type                = input_type
  , typename coalesced_output_type      = output_type
>
struct apply_wrapper
{
    cuda::function<void (coalesced_input_type const*, coalesced_output_type*, unsigned int)> apply;
    static apply_wrapper const kernel;
};

} // namespace algorithm
} // namespace gpu
} // namespace halmd

#endif /* ! HALMD_ALGORITHM_GPU_APPLY_KERNEL_HPP */
