/* cuda_wrapper/mem.hpp
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

#ifndef CUDA_DRIVER_MEM_HPP
#define CUDA_DRIVER_MEM_HPP

#include <boost/shared_ptr.hpp>
#include <cuda.h>
#include <cuda_wrapper/driver/error.hpp>
#include <string>

namespace cuda { namespace driver { namespace mem
{

/**
 * get total memory in bytes of current CUDA context
 */
inline unsigned int total()
{
    unsigned int free = 0, total = 0;
    CU_CALL(cuMemGetInfo(&free, &total));
    return total;
}

/**
 * get allocated memory in bytes of current CUDA context
 */
inline unsigned int used()
{
    unsigned int free = 0, total = 0;
    CU_CALL(cuMemGetInfo(&free, &total));
    return (total - free);
}

/**
 * get available memory in bytes of current CUDA context
 */
inline unsigned int free()
{
    unsigned int free = 0, total = 0;
    CU_CALL(cuMemGetInfo(&free, &total));
    return free;
}

}}} // namespace cuda::driver::mem

#endif /* ! CUDA_DRIVER_MEM_HPP */
