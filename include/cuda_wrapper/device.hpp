/* cuda_wrapper/device.hpp
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

#ifndef CUDA_DEVICE_HPP
#define CUDA_DEVICE_HPP

#include <cuda/cuda_runtime.h>
#include <cuda/cuda.h>
#include <cuda_wrapper/error.hpp>
#include <string>

namespace cuda
{

/**
 * CUDA device management
 */
class device
{
public:
    /**
     * returns number of devices available for execution
     */
    static int count()
    {
	int count;
	CUDA_CALL(cudaGetDeviceCount(&count));
	return count;
    }

    /**
     * set device on which the active host thread executes device code
     */
    static void set(int dev)
    {
	CUDA_CALL(cudaSetDevice(dev));
    }

    /**
     * get device on which the active host thread executes device code
     */
    static int get()
    {
	int dev;
	CUDA_CALL(cudaGetDevice(&dev));
	return dev;
    }

    /**
     * CUDA device properties
     */
    class properties
    {
    private:
	cudaDeviceProp prop;

    public:
	/**
	 * empty initializer
	 */
	properties() {}

	/**
	 * retrieve properties of given device
	 */
	properties(int dev)
	{
	    CUDA_CALL(cudaGetDeviceProperties(&prop, dev));
	}

	/**
	 * ASCII string identifying the device
	 */
	std::string name()
	{
	    return prop.name;
	}

	/**
	 * total amount of global memory available on the device in bytes
	 */
	size_t total_global_mem()
	{
	    return prop.totalGlobalMem;
	}

	/**
	 * total amount of shared memory available per block in bytes
	 */
	size_t shared_mem_per_block()
	{
	    return prop.sharedMemPerBlock;
	}

	/**
	 * total number of registers available per block
	 */
	size_t regs_per_block()
	{
	    return prop.regsPerBlock;
	}

	/**
	 * wrap size
	 */
	size_t warp_size()
	{
	    return prop.warpSize;
	}

	/**
	 * maximum allowed memory allocation pitch
	 */
	size_t mem_pitch()
	{
	    return prop.memPitch;
	}

	/**
	 * maximum number of threads per block
	 */
	unsigned int max_threads_per_block()
	{
	    return prop.maxThreadsPerBlock;
	}

	/**
	 * maximum sizes of each dimension of a block
	 */
	dim3 max_threads_dim()
	{
	    return dim3(prop.maxThreadsDim[0], prop.maxThreadsDim[1], prop.maxThreadsDim[2]);
	}

	/**
	 * maximum sizes of each dimension of a grid
	 */
	dim3 max_grid_size()
	{
	    return dim3(prop.maxGridSize[0], prop.maxGridSize[1], prop.maxGridSize[2]);
	}

	/**
	 * total amount of constant memory available on the device in bytes
	 */
	size_t total_const_mem()
	{
	    return prop.totalConstMem;
	}

	/**
	 * major revision number of device's compute capatibility
	 */
	unsigned int major()
	{
	    return prop.major;
	}

	/**
	 * minor revision number of device's compute capatibility
	 */
	unsigned int minor()
	{
	    return prop.minor;
	}

	/**
	 * clock frequency in kHz
	 */
	unsigned int clock_rate()
	{
	    return prop.clockRate;
	}

	/**
	 * texture alignment requirement
	 */
	size_t texture_alignment()
	{
	    return prop.textureAlignment;
	}
    };

    /**
     * get total memory in bytes for given device
     */
    static unsigned int mem_get_total(int dev)
    {
	unsigned int free = 0, total = 0;
	_mem_get_info(&free, &total, dev);
	return total;
    }

    /**
     * get allocated memory in bytes for given device
     */
    static unsigned int mem_get_used(int dev)
    {
	unsigned int free = 0, total = 0;
	_mem_get_info(&free, &total, dev);
	return (total - free);
    }

    /**
     * get available memory in bytes for given device
     */
    static unsigned int mem_get_free(int dev)
    {
	unsigned int free = 0, total = 0;
	_mem_get_info(&free, &total, dev);
	return free;
    }

private:
    /**
     * get free and total memory in the current context
     */
    static void _mem_get_info(unsigned int* free, unsigned int* total, int dev)
    {
	CUcontext cuctx;
	CUdevice cudev;

	/* create CUDA context for device */
	CU_CALL(cuInit(0));
	CU_CALL(cuDeviceGet(&cudev, dev));
	CU_CALL(cuCtxCreate(&cuctx, 0, cudev));
	/* query memory info */
	CU_CALL(cuMemGetInfo(free, total));
	/* restore previous context, if any */
	CU_CALL(cuCtxPopCurrent(NULL));
	/* destroy CUDA context */
	CU_CALL(cuCtxDestroy(cuctx));
    }
};

}

#endif /* ! CUDA_DEVICE_HPP */
