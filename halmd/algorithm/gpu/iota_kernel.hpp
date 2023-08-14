/*
 * Copyright © 2012 Peter Colberg
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

#ifndef HALMD_ALGORITHM_GPU_IOTA_KERNEL_HPP
#define HALMD_ALGORITHM_GPU_IOTA_KERNEL_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>

namespace halmd {
namespace detail {

class iota_kernel
{
public:
    cuda::function<void (unsigned int*, unsigned int, unsigned int)> iota;

    static iota_kernel& get()
    {
        return kernel;
    }

private:
    static iota_kernel kernel;
};

} // namespace detail
} // namespace halmd

#endif /* ! HALMD_ALGORITHM_GPU_IOTA_KERNEL_HPP */
