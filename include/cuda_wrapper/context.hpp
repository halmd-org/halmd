/* cuda_wrapper/context.hpp
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

#ifndef CUDA_CONTEXT_HPP
#define CUDA_CONTEXT_HPP

#include <boost/shared_ptr.hpp>
#include <cuda/cuda.h>
#include <cuda_wrapper/error.hpp>
#include <string>

namespace cuda
{

/**
 * CUDA context management
 */
class context
{
public:
    /**
     * create a CUDA context and associate it with the calling thread
     */
    context(int ordinal, unsigned int flags = 0) : m_ctx(new CUcontext)
    {
	CUdevice dev;
	CU_CALL(cuInit(0));
	CU_CALL(cuDeviceGet(&dev, ordinal));
	CU_CALL(cuCtxCreate(&(*m_ctx), flags, dev));
    }

    /**
     * detaches a CUDA context, and destroys it if usage count is zero
     */
    ~context()
    {
	// discard return value to guarantee no-throw deconstructor
	cuCtxDetach(*m_ctx);
    }

    /**
     * increment usage count of a context
     */
    context(context const& ctx) : m_ctx(ctx.m_ctx)
    {
	CU_CALL(cuCtxAttach(&(*m_ctx), 0));
    }

    context& operator=(context const& ctx)
    {
	CU_CALL(cuCtxDetach(*m_ctx));
	m_ctx = ctx.m_ctx;
	CU_CALL(cuCtxAttach(&(*m_ctx), 0));
	return *this;
    }

    /**
     * returns ordinal of the current context's device
     */
    static CUdevice device()
    {
	CUdevice dev;
	CU_CALL(cuCtxGetDevice(&dev));
	return dev;
    }

    /**
     * block for a context's tasks to complete
     */
    static void synchronize()
    {
	CU_CALL(cuCtxSynchronize());
    }

private:
    boost::shared_ptr<CUcontext> m_ctx;
};

} // namespace cuda

#endif /* ! CUDA_CONTEXT_HPP */
