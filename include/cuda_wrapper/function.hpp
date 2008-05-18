/* CUDA device function execution
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

#ifndef BOOST_PP_IS_ITERATING

    #ifndef CUDA_FUNCTION_HPP
    #define CUDA_FUNCTION_HPP

    #include <boost/preprocessor/iteration/iterate.hpp>
    #include <boost/preprocessor/repetition/enum_params.hpp>
    #include <boost/preprocessor/repetition/enum_binary_params.hpp>
    #include <boost/preprocessor/repetition/repeat.hpp>
    #include <cuda/cuda_runtime.h>
    #ifndef __CUDACC__
    #include <cuda_wrapper/error.hpp>
    #include <cuda_wrapper/stream.hpp>
    #endif


    /* maximum number of arguments passed to device functions */
    #ifndef CUDA_FUNCTION_MAX_ARGS
    #define CUDA_FUNCTION_MAX_ARGS 10
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
	dim3 grid;
	/* block dimensions */
	dim3 block;
	/* FIXME store useful numbers (no. of threads per grid/block) */

	config()
	{
	}

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

    } // namespace cuda


    #define BOOST_PP_FILENAME_1 <cuda_wrapper/function.hpp>
    #define BOOST_PP_ITERATION_LIMITS (1, CUDA_FUNCTION_MAX_ARGS)
    #include BOOST_PP_ITERATE()

    #endif /* ! CUDA_FUNCTION_HPP */

#else /* ! BOOST_PP_IS_ITERATING */

    namespace cuda
    {

    template <typename T>
    class function;

    /**
     * CUDA kernel execution wrapper for n-ary device function
     */
    template <BOOST_PP_ENUM_PARAMS(BOOST_PP_ITERATION(), typename T)>
    class function<void (BOOST_PP_ENUM_PARAMS(BOOST_PP_ITERATION(), T))>
    {
    public:
	typedef void T (BOOST_PP_ENUM_PARAMS(BOOST_PP_ITERATION(), T));

    public:
	function(T *entry) : entry(entry) {}

    #ifndef __CUDACC__

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

	/**
	 * execute kernel
	 */
	void operator()(BOOST_PP_ENUM_BINARY_PARAMS(BOOST_PP_ITERATION(), const T, &x))
	{
	    size_t offset = 0;
    #define SETUP_ARGUMENT(z, n, x) setup_argument(x##n, &offset);
	    BOOST_PP_REPEAT(BOOST_PP_ITERATION(), SETUP_ARGUMENT, x)
    #undef SETUP_ARGUMENT
	    CUDA_CALL(cudaLaunch(reinterpret_cast<const char *>(entry)));
	}

    private:
	/**
	 * push arbitrary argument into argument passing area
	 */
	template <typename U>
	static void setup_argument(const U& arg, size_t *offset)
	{
	    /* respect alignment requirements of passed argument */
	    if (0 != *offset % __alignof(U)) {
		*offset += __alignof(U) - *offset % __alignof(U);
	    }

	    CUDA_CALL(cudaSetupArgument(&arg, sizeof(U), *offset));

	    /* advance argument offset for next call */
	    *offset += sizeof(U);
	}

    #endif /* ! __CUDACC__ */

    private:
	T *entry;
    };

    } // namespace cuda

#endif /* ! BOOST_PP_IS_ITERATING */
