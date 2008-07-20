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
    #define CUDA_FUNCTION_MAX_ARGS 20
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

    #ifndef __CUDACC__

    /**
     * configure execution parameters
     */
    __inline__ void configure(dim3 const& grid, dim3 const& block, size_t shared_mem = 0)
    {
	CUDA_CALL(cudaConfigureCall(grid, block, shared_mem, 0));
    }

    #ifdef CUDA_WRAPPER_ASYNC_API
    /**
     * configure execution parameters
     */
    __inline__ void configure(dim3 const& grid, dim3 const& block, stream& stream)
    {
	CUDA_CALL(cudaConfigureCall(grid, block, 0, stream.data()));
    }

    /**
     * configure execution parameters
     */
    __inline__ void configure(dim3 const& grid, dim3 const& block, size_t shared_mem, stream& stream)
    {
	CUDA_CALL(cudaConfigureCall(grid, block, shared_mem, stream.data()));
    }
    #endif /* CUDA_WRAPPER_ASYNC_API */

    /**
     * CUDA device function argument wrapper
     */
    template <typename T>
    class arg
    {
    public:
	arg(T const& a) : a(a) {}

	/**
	 * returns const reference to arbitrary argument
	 */
	T const& operator*()
	{
	    return a;
	}

    private:
	T const& a;
    };

    template <typename T>
    class arg<T*>
    {
    public:
	arg(vector<T>& a) : a(a) {}

	/**
	 * returns pointer to CUDA device memory array
	 */
	T* operator*()
	{
	    return a.data();
	}

    private:
	vector<T>& a;
    };

    template <typename T>
    class arg<T const*>
    {
    public:
	arg(vector<T> const& a) : a(a) {}

	/**
	 * returns const pointer to CUDA device memory array
	 */
	T const* operator*()
	{
	    return a.data();
	}

    private:
	vector<T> const& a;
    };

    #endif /* ! __CUDACC__ */

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
	 * execute kernel
	 */
	void operator()(BOOST_PP_ENUM_BINARY_PARAMS(BOOST_PP_ITERATION(), arg<T, > x))
	{
	    // properly align CUDA device function arguments
	    struct {
		#define DECL_ARG(z, n, x) T##n x##n;
		BOOST_PP_REPEAT(BOOST_PP_ITERATION(), DECL_ARG, a)
		#undef DECL_ARG
	    } args = {
		BOOST_PP_ENUM_PARAMS(BOOST_PP_ITERATION(), *x)
	    };
	    // push aligned arguments onto CUDA execution stack
	    CUDA_CALL(cudaSetupArgument(&args, sizeof(args), 0));
	    // launch CUDA device function
	    CUDA_CALL(cudaLaunch(reinterpret_cast<const char *>(entry)));
	}

    #endif /* ! __CUDACC__ */

    private:
	T *entry;
    };

    } // namespace cuda

#endif /* ! BOOST_PP_IS_ITERATING */
