/* cuda_wrapper/thread.hpp
 *
 * Copyright (C) 2007  Peter Colberg
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

/*
 * CUDA thread management
 */

#ifndef CUDA_THREAD_HPP
#define CUDA_THREAD_HPP

#include <cuda/cuda_runtime.h>
#include <cuda_wrapper/error.hpp>

namespace cuda
{

/*
 * blocks until the device has completed all preceding requested tasks
 */
static void thread_synchronize()
{
    CUDA_CALL(cudaThreadSynchronize());
}

/*
 * cleans up all runtime-related resources associated with calling thread
 */
static void thread_exit()
{
    CUDA_CALL(cudaThreadExit());
}

}

#endif /* ! CUDA_THREAD_HPP */
