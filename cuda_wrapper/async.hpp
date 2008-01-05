/* cuda_wrapper/async.hpp
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
 * Asynchronous CUDA runtime API wrapper classes
 */

#ifndef CUDA_ASYNC_HPP
#define CUDA_ASYNC_HPP

/* requires CUDA runtime version >= 1.1 */
#if (CUDART_VERSION >= 1010)
#define CUDA_ASYNC_API
#endif

#ifdef CUDA_ASYNC_API

#include <cuda_runtime.h>
#include <cuda_wrapper/error.hpp>

/*
 * CUDA stream wrapper class
 */
class cuda_stream
{
    friend class cuda_event;

protected:
    cudaStream_t _stream;

public:
    /*
     * creates a stream
     */
    cuda_stream()
    {
	CUDA_CALL(cudaStreamCreate(&_stream));
    }

    /*
     * destroys the stream
     */
    ~cuda_stream()
    {
	CUDA_CALL(cudaStreamDestroy(_stream));
    }

    /*
     * blocks until the device has completed all operations in the stream
     */
    void synchronize()
    {
	CUDA_CALL(cudaStreamSynchronize(_stream));
    }

    /*
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
    cuda_stream(const cuda_stream&);
    // disable default assignment operator
    cuda_stream& operator=(const cuda_stream&);
};


/*
 * CUDA event wrapper class
 */
class cuda_event
{
protected:
    cudaEvent_t _event;

public:
    /*
     * creates an event
     */
    cuda_event()
    {
	CUDA_CALL(cudaEventCreate(&_event));
    }

    /*
     * destroys the event
     */
    ~cuda_event()
    {
	CUDA_CALL(cudaEventDestroy(_event));
    }

    /*
     * records an event
     *
     * after all preceding operations in the CUDA context have been completed
     */
    void record()
    {
	CUDA_CALL(cudaEventRecord(_event, 0));
    }

    /*
     * records an event
     *
     * after all preceding operations in the stream have been completed
     */
    void record(const cuda_stream& stream)
    {
	CUDA_CALL(cudaEventRecord(_event, stream._stream));
    }

    /*
     * blocks until the event has actually been recorded
     */
    void synchronize()
    {
	CUDA_CALL(cudaEventSynchronize(_event));
    }

    /*
     * checks if the event has actually been recorded
     *
     * WARNING: this function will not detect kernel launch failures
     */
    bool query()
    {
	cudaError_t err = cudaEventQuery(_event);
	if (cudaSuccess == err)
	    return true;
	else if (cudaErrorNotReady == err)
	    return false;
	CUDA_ERROR(err);
    }

    /*
     * computes the elapsed time between two events
     *
     * (in milliseconds with a resolution of around 0.5 microseconds)
     */
    float operator-(const cuda_event &start)
    {
	float time;
	CUDA_CALL(cudaEventElapsedTime(&time, start._event, _event));
	return time;
    }

private:
    // disable default copy constructor
    cuda_event(const cuda_event&);
    // disable default assignment operator
    cuda_event& operator=(const cuda_event&);
};

#endif /* CUDA_ASYNC_API */

#endif /* ! CUDA_ASYNC_HPP */
