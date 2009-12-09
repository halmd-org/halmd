/* cuda_wrapper/context.hpp
 *
 * Copyright (C) 2009  Peter Colberg
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

#ifndef CUDA_DRIVER_CONTEXT_HPP
#define CUDA_DRIVER_CONTEXT_HPP

#include <boost/shared_ptr.hpp>
#include <cuda.h>
#include <cuda_wrapper/driver/error.hpp>
#include <string>

namespace cuda { namespace driver
{

/**
 * CUDA context management
 */
class context
{
public:
    /**
     * transform current context into floating context
     */
    class floating : boost::noncopyable
    {
    public:
        /**
         * pop current context from CUDA context stack
         */
        floating()
        {
            // attach to current context
            CU_CALL(cuCtxAttach(&m_ctx, 0));
            // copy context pointer
            CUcontext ctx(m_ctx);
            // set internal usage count
            m_usage = 0;

            while (ctx == m_ctx) {
                // decrement usage count of current context
                CU_CALL(cuCtxDetach(ctx));
                // requires context usage count of 1
                CU_CALL(cuCtxPopCurrent(NULL));
                try {
                    // is there a current context?
                    CU_CALL(cuCtxAttach(&ctx, 0));
                }
                catch (cuda::driver::error const&) {
                    // no current context
                    break;
                }
                // decrement usage count of current context
                CU_CALL(cuCtxDetach(ctx));
                // increment internal usage count
                m_usage++;
            }
        }

        /**
         * push floating context onto CUDA context stack
         */
        ~floating()
        {
            CU_CALL(cuCtxPushCurrent(m_ctx));
            // restore usage count of context
            while (m_usage) {
                CU_CALL(cuCtxAttach(&m_ctx, 0));
                m_usage--;
            }
        }

    private:
        CUcontext m_ctx;
        size_t m_usage;
    };

public:
    /**
     * create a CUDA context and associate it with the calling thread
     */
#if (CUDART_VERSION >= 2020)
    context(int ordinal, unsigned int flags = CU_CTX_MAP_HOST) : m_ctx(new CUcontext)
#else
    context(int ordinal, unsigned int flags = 0) : m_ctx(new CUcontext)
#endif
    {
        CUdevice dev;
        CU_CALL(cuInit(0));
        CU_CALL(cuDeviceGet(&dev, ordinal));
        CU_CALL(cuCtxCreate(&(*m_ctx), flags, dev));
    }

    /**
     * attaches to the current CUDA context
     */
    context() : m_ctx(new CUcontext)
    {
        CU_CALL(cuCtxAttach(&(*m_ctx), 0));
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

}} // namespace cuda::driver

#endif /* ! CUDA_DRIVER_CONTEXT_HPP */
