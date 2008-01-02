/* cuda_wrapper/error.hpp
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
 * CUDA runtime error checking
 */

#ifndef CUDA_ERROR_HPP
#define CUDA_ERROR_HPP

#include <cuda/cuda_runtime.h>


#define CUDA_ERROR(err) throw cuda::error(err)

#define CUDA_CALL(x)							\
    do {								\
	cudaError_t err;						\
	if (cudaSuccess != (err = x)) {					\
	    CUDA_ERROR(err);						\
	}								\
    } while(0)


namespace cuda
{

/*
 * CUDA error handling
 */
class error
{
public:
    /* CUDA error */
    const cudaError_t err;

    error(cudaError_t err): err(err)
    {
    }

    /*
     * returns a message string for the CUDA error
     */
    const char* what() const throw()
    {
	return cudaGetErrorString(err);
    }
};

}

#endif /* ! CUDA_ERROR_HPP */
