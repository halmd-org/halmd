/* cuda_wrapper/driver/version.hpp
 *
 * Copyright (C) 2009  Peter Colberg
 *
 * This program is free software: you can redistribute it and/or modify
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

#ifndef CUDA_DRIVER_VERSION_HPP
#define CUDA_DRIVER_VERSION_HPP

#include <cuda.h>
#include <cuda_wrapper/driver/error.hpp>

namespace cuda { namespace driver
{

#if (CUDA_VERSION >= 2020)

/**
 * returns CUDA driver library version
 */
inline int version()
{
    int v;
    CU_CALL(cuDriverGetVersion(&v));
    return v;
}

#endif /* CUDA_VERSION >= 2020 */

}} // namespace cuda::driver

#endif /* ! CUDA_DRIVER_VERSION_HPP */
