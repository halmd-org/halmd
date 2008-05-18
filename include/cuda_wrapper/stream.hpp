/* cuda_wrapper/stream.hpp
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

#ifndef CUDA_STREAM_HPP
#define CUDA_STREAM_HPP

#include <cuda/cuda_runtime.h>
#include <cuda_wrapper/error.hpp>


namespace cuda
{

template <typename T>
class vector;

namespace host
{

template <typename T, typename Alloc>
class vector;

}

#ifdef CUDA_WRAPPER_ASYNC_API

/**
 * CUDA stream wrapper class
 */
class stream
{
    friend class event;

    template <typename T>
    friend class vector;
    template <typename T, typename Alloc>
    friend class host::vector;

    template <typename T>
    friend class function;

protected:
    cudaStream_t _stream;

public:
    /**
     * creates a stream
     */
    stream()
    {
	CUDA_CALL(cudaStreamCreate(&_stream));
    }

    /**
     * destroys the stream
     */
    ~stream()
    {
	CUDA_CALL(cudaStreamDestroy(_stream));
    }

    /**
     * blocks until the device has completed all operations in the stream
     */
    void synchronize()
    {
	CUDA_CALL(cudaStreamSynchronize(_stream));
    }

    /**
     * checks if the device has completed all operations in the stream
     *
     * WARNING: this function will not detect kernel launch failures
     */
    bool query()
    {
	cudaError_t err = cudaStreamQuery(_stream);
	if (cudaSuccess == err)
	    return true;
	else if (cudaErrorNotReady == err)
	    return false;
	CUDA_ERROR(err);
    }

private:
    // disable default copy constructor
    stream(const stream&);
    // disable default assignment operator
    stream& operator=(const stream&);
};

#endif /* CUDA_WRAPPER_ASYNC_API */

}

#endif /* ! CUDA_STREAM_HPP */
