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

#include <cuda_runtime.h>
#include <cuda.h>
#include <exception>


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
class error : public std::exception
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

#define CU_ERROR(err) throw cu::error(err)

#define CU_CALL(x)							\
    do {								\
	CUresult err;							\
	if (CUDA_SUCCESS != (err = x)) {					\
	    CU_ERROR(err);						\
	}								\
    } while(0)


namespace cu
{

class error : public std::exception
{
public:
    const CUresult err;

    error(CUresult err): err(err)
    {
    }

    /*
     * returns a message string for the CUDA error
     */
    const char* what() const throw()
    {
	switch (err) {
	  case CUDA_SUCCESS:
	    return "No errors";
	  case CUDA_ERROR_INVALID_VALUE:
	    return "Invalid value";
	  case CUDA_ERROR_OUT_OF_MEMORY:
	    return "Out of memory";
	  case CUDA_ERROR_NOT_INITIALIZED:
	    return "Driver not initialized";
	  case CUDA_ERROR_DEINITIALIZED:
	    return "Driver deinitialized";
	  case CUDA_ERROR_NO_DEVICE:
	    return "No CUDA-capable device available";
	  case CUDA_ERROR_INVALID_DEVICE:
	    return "Invalid device";
	  case CUDA_ERROR_INVALID_IMAGE:
	    return "Invalid kernel image";
	  case CUDA_ERROR_INVALID_CONTEXT:
	    return "Invalid context";
	  case CUDA_ERROR_CONTEXT_ALREADY_CURRENT:
	    return "Context already current";
	  case CUDA_ERROR_MAP_FAILED:
	    return "Map failed";
	  case CUDA_ERROR_UNMAP_FAILED:
	    return "Unmap failed";
	  case CUDA_ERROR_ARRAY_IS_MAPPED:
	    return "Array is mapped";
	  case CUDA_ERROR_ALREADY_MAPPED:
	    return "Already mapped";
	  case CUDA_ERROR_NO_BINARY_FOR_GPU:
	    return "No binary for GPU";
	  case CUDA_ERROR_ALREADY_ACQUIRED:
	    return "Already acquired";
	  case CUDA_ERROR_NOT_MAPPED:
	    return "Not mapped";
	  case CUDA_ERROR_INVALID_SOURCE:
	    return "Invalid source";
	  case CUDA_ERROR_FILE_NOT_FOUND:
	    return "File not found";
	  case CUDA_ERROR_INVALID_HANDLE:
	    return "Invalid handle";
	  case CUDA_ERROR_NOT_FOUND:
	    return "Not found";
	  case CUDA_ERROR_NOT_READY:
	    return "CUDA not ready";
	  case CUDA_ERROR_LAUNCH_FAILED:
	    return "Launch failed";
	  case CUDA_ERROR_LAUNCH_OUT_OF_RESOURCES:
	    return "Launch exceeded resources";
	  case CUDA_ERROR_LAUNCH_TIMEOUT:
	    return "Launch exceeded timeout";
	  case CUDA_ERROR_LAUNCH_INCOMPATIBLE_TEXTURING:
	    return "Launch with incompatible texturing";
	  default:
	    return "Unknown error";
	}
    }
};

}

#endif /* ! CUDA_ERROR_HPP */
