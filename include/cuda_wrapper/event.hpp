/* cuda_wrapper/event.hpp
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

#ifndef CUDA_EVENT_HPP
#define CUDA_EVENT_HPP

#include <cuda_runtime.h>
#include <cuda_wrapper/error.hpp>
#include <cuda_wrapper/stream.hpp>

namespace cuda
{

#ifdef CUDA_WRAPPER_ASYNC_API

/**
 * CUDA event wrapper class
 */
class event
{
protected:
    cudaEvent_t _event;

public:
    /**
     * creates an event
     */
    event()
    {
	CUDA_CALL(cudaEventCreate(&_event));
    }

    /**
     * destroys the event
     */
    ~event()
    {
	CUDA_CALL(cudaEventDestroy(_event));
    }

    /**
     * records an event
     *
     * after all preceding operations in the CUDA context have been completed
     */
    void record()
    {
	CUDA_CALL(cudaEventRecord(_event, 0));
    }

    /**
     * records an event
     *
     * after all preceding operations in the stream have been completed
     */
    void record(const stream& stream)
    {
	CUDA_CALL(cudaEventRecord(_event, stream._stream));
    }

    /**
     * blocks until the event has actually been recorded
     */
    void synchronize()
    {
	CUDA_CALL(cudaEventSynchronize(_event));
    }

    /**
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

    /**
     * computes the elapsed time between two events
     *
     * (in milliseconds with a resolution of around 0.5 microseconds)
     */
    float operator-(const event &start)
    {
	float time;
	CUDA_CALL(cudaEventElapsedTime(&time, start._event, _event));
	return time;
    }

private:
    // disable default copy constructor
    event(const event&);
    // disable default assignment operator
    event& operator=(const event&);
};

#endif /* CUDA_WRAPPER_ASYNC_API */

}

#endif /* ! CUDA_EVENT_HPP */
