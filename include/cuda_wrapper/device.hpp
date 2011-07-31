/* cuda_wrapper/device.hpp
 *
 * Copyright (C) 2007  Peter Colberg
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

#ifndef CUDA_DEVICE_HPP
#define CUDA_DEVICE_HPP

#include <cuda_runtime.h>
#include <string>

#include <cuda_wrapper/error.hpp>

namespace cuda {

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
        std::string name() const
        {
            return prop.name;
        }

        /**
         * total amount of global memory available on the device in bytes
         */
        size_t total_global_mem() const
        {
            return prop.totalGlobalMem;
        }

        /**
         * total amount of shared memory available per block in bytes
         */
        size_t shared_mem_per_block() const
        {
            return prop.sharedMemPerBlock;
        }

        /**
         * total number of registers available per block
         */
        size_t regs_per_block() const
        {
            return prop.regsPerBlock;
        }

        /**
         * wrap size
         */
        size_t warp_size() const
        {
            return prop.warpSize;
        }

        /**
         * maximum allowed memory allocation pitch
         */
        size_t mem_pitch() const
        {
            return prop.memPitch;
        }

        /**
         * maximum number of threads per block
         */
        unsigned int max_threads_per_block() const
        {
            return prop.maxThreadsPerBlock;
        }

        /**
         * maximum sizes of each dimension of a block
         */
        dim3 max_threads_dim() const
        {
            return dim3(prop.maxThreadsDim[0], prop.maxThreadsDim[1], prop.maxThreadsDim[2]);
        }

        /**
         * maximum sizes of each dimension of a grid
         */
        dim3 max_grid_size() const
        {
            return dim3(prop.maxGridSize[0], prop.maxGridSize[1], prop.maxGridSize[2]);
        }

        /**
         * total amount of constant memory available on the device in bytes
         */
        size_t total_const_mem() const
        {
            return prop.totalConstMem;
        }

        /**
         * major revision number of device's compute capatibility
         */
        unsigned int major() const
        {
            return prop.major;
        }

        /**
         * minor revision number of device's compute capatibility
         */
        unsigned int minor() const
        {
            return prop.minor;
        }

        /**
         * clock frequency in kHz
         */
        unsigned int clock_rate() const
        {
            return prop.clockRate;
        }

        /**
         * texture alignment requirement
         */
        size_t texture_alignment() const
        {
            return prop.textureAlignment;
        }

#if (CUDART_VERSION >= 2000)
        /**
         * asynchronous kernel and memory operations capability
         */
        int device_overlap() const
        {
            return prop.deviceOverlap;
        }

        /**
         * number of multiprocessors
         */
        int multi_processor_count() const
        {
            return prop.multiProcessorCount;
        }
#endif /* CUDART_VERSION >= 2000 */
    };
};

} // namespace cuda

#endif /* ! CUDA_DEVICE_HPP */
