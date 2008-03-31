/* cuda_wrapper/function.hpp
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

#ifndef CUDA_FUNCTION_HPP
#define CUDA_FUNCTION_HPP

#include <cuda_runtime.h>
#ifndef __CUDACC__
#include <cuda_wrapper/error.hpp>
#include <cuda_wrapper/stream.hpp>
#endif


namespace cuda
{

/*
 * CUDA execution configuration
 */
class config
{
public:
    /* grid dimensions */
    const dim3 grid;
    /* block dimensions */
    const dim3 block;
    /* FIXME store useful numbers (no. of threads per grid/block) */

    config(dim3 grid, dim3 block) : grid(grid), block(block)
    {
	/* FIXME store useful numbers (no. of threads per grid/block) */
    }

    size_t threads() const
    {
	return grid.y * grid.x * block.z * block.y * block.x;
    }

    size_t blocks_per_grid() const
    {
	return grid.y * grid.x;
    }

    size_t threads_per_block() const
    {
	return block.z * block.y * block.x;
    }
};


/**
 * CUDA kernel execution wrapper base class
 */
class _function_base
{
#ifndef __CUDACC__
public:
    /**
     * configure execution parameters
     */
    static void configure(const config& dim, size_t shared_mem = 0)
    {
	CUDA_CALL(cudaConfigureCall(dim.grid, dim.block, shared_mem, 0));
    }

#ifdef CUDA_WRAPPER_ASYNC_API

    /**
     * configure execution parameters
     */
    static void configure(const config& dim, stream& stream)
    {
	CUDA_CALL(cudaConfigureCall(dim.grid, dim.block, 0, stream._stream));
    }

    /**
     * configure execution parameters
     */
    static void configure(const config& dim, size_t shared_mem, stream& stream)
    {
	CUDA_CALL(cudaConfigureCall(dim.grid, dim.block, shared_mem, stream._stream));
    }

#endif /* CUDA_WRAPPER_ASYNC_API */

protected:
    /**
     * push arbitrary argument into argument passing area
     */
    template <typename T>
    static void setup_argument(const T& arg, size_t *offset)
    {
	/* respect alignment requirements of passed argument */
	if (0 != *offset % __alignof(T)) {
	    *offset += __alignof(T) - *offset % __alignof(T);
	}

	CUDA_CALL(cudaSetupArgument(&arg, sizeof(T), *offset));

	/* advance argument offset for next call */
	*offset += sizeof(T);
    }

    /**
     * launch kernel
     */
    template <typename T>
    static void launch(T *entry)
    {
	CUDA_CALL(cudaLaunch(reinterpret_cast<const char *>(entry)));
    }
#endif /* ! __CUDACC__ */
};


/**
 * CUDA kernel execution wrapper
 */
template <typename T>
class function;


/**
 * CUDA kernel execution wrapper for unary device function
 */
template <typename T0>
class function<void (T0)> : public _function_base
{
protected:
    typedef void T (T0);
    T *entry;

public:
    function(T *entry) : entry(entry) {}

#ifndef __CUDACC__
    void operator()(T0 x0)
    {
	size_t offset = 0;
	setup_argument(x0, &offset);
	launch(entry);
    }
#endif /* ! __CUDACC__ */
};


/**
 * CUDA kernel execution wrapper for binary device function
 */
template <typename T0, typename T1>
class function<void (T0, T1)> : public _function_base
{
protected:
    typedef void T (T0, T1);
    T *entry;

public:
    function(T *entry) : entry(entry) {}

#ifndef __CUDACC__
    void operator()(T0 x0, T1 x1)
    {
	size_t offset = 0;
	setup_argument(x0, &offset);
	setup_argument(x1, &offset);
	launch(entry);
    }
#endif /* ! __CUDACC__ */
};


/**
 * CUDA kernel execution wrapper for ternary device function
 */
template <typename T0, typename T1, typename T2>
class function<void (T0, T1, T2)> : public _function_base
{
protected:
    typedef void T (T0, T1, T2);
    T *entry;

public:
    function(T *entry) : entry(entry) {}

#ifndef __CUDACC__
    void operator()(T0 x0, T1 x1, T2 x2)
    {
	size_t offset = 0;
	setup_argument(x0, &offset);
	setup_argument(x1, &offset);
	setup_argument(x2, &offset);
	launch(entry);
    }
#endif /* ! __CUDACC__ */
};


/**
 * CUDA kernel execution wrapper for quaternary device function
 */
template <typename T0, typename T1, typename T2, typename T3>
class function<void (T0, T1, T2, T3)> : public _function_base
{
protected:
    typedef void T (T0, T1, T2, T3);
    T *entry;

public:
    function(T *entry) : entry(entry) {}

#ifndef __CUDACC__
    void operator()(T0 x0, T1 x1, T2 x2, T3 x3)
    {
	size_t offset = 0;
	setup_argument(x0, &offset);
	setup_argument(x1, &offset);
	setup_argument(x2, &offset);
	setup_argument(x3, &offset);
	launch(entry);
    }
#endif /* ! __CUDACC__ */
};


/**
 * CUDA kernel execution wrapper for quinary device function
 */
template <typename T0, typename T1, typename T2, typename T3, typename T4>
class function<void (T0, T1, T2, T3, T4)> : public _function_base
{
protected:
    typedef void T (T0, T1, T2, T3, T4);
    T *entry;

public:
    function(T *entry) : entry(entry) {}

#ifndef __CUDACC__
    void operator()(T0 x0, T1 x1, T2 x2, T3 x3, T4 x4)
    {
	size_t offset = 0;
	setup_argument(x0, &offset);
	setup_argument(x1, &offset);
	setup_argument(x2, &offset);
	setup_argument(x3, &offset);
	setup_argument(x4, &offset);
	launch(entry);
    }
#endif /* ! __CUDACC__ */
};


/**
 * CUDA kernel execution wrapper for senary device function
 */
template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5>
class function<void (T0, T1, T2, T3, T4, T5)> : public _function_base
{
protected:
    typedef void T (T0, T1, T2, T3, T4, T5);
    T *entry;

public:
    function(T *entry) : entry(entry) {}

#ifndef __CUDACC__
    void operator()(T0 x0, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5)
    {
	size_t offset = 0;
	setup_argument(x0, &offset);
	setup_argument(x1, &offset);
	setup_argument(x2, &offset);
	setup_argument(x3, &offset);
	setup_argument(x4, &offset);
	setup_argument(x5, &offset);
	launch(entry);
    }
#endif /* ! __CUDACC__ */
};


/**
 * CUDA kernel execution wrapper for septenary device function
 */
template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
class function<void (T0, T1, T2, T3, T4, T5, T6)> : public _function_base
{
protected:
    typedef void T (T0, T1, T2, T3, T4, T5, T6);
    T *entry;

public:
    function(T *entry) : entry(entry) {}

#ifndef __CUDACC__
    void operator()(T0 x0, T1 x1, T2 x2, T3 x3, T4 x4, T5 x5, T6 x6)
    {
	size_t offset = 0;
	setup_argument(x0, &offset);
	setup_argument(x1, &offset);
	setup_argument(x2, &offset);
	setup_argument(x3, &offset);
	setup_argument(x4, &offset);
	setup_argument(x5, &offset);
	setup_argument(x6, &offset);
	launch(entry);
    }
#endif /* ! __CUDACC__ */
};

}

#endif /* ! CUDA_FUNCTION_HPP */
