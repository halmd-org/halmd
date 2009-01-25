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

#if !BOOST_PP_IS_ITERATING

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

    /* enable function overloading for varying number of arguments */
    #ifndef CUDA_FUNCTION_VARY_ARGS
    #define CUDA_FUNCTION_VARY_ARGS 0	/* set to 1 for extended coffee break… */
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

    #if (CUDART_VERSION >= 1010)

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

    #endif /* CUDART_VERSION >= 1010 */

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
	arg(vector<T>& a) : a(a.data()) {}
	arg(T* a) : a(a) {}

	/**
	 * returns pointer to CUDA device memory array
	 */
	T* operator*()
	{
	    return a;
	}

    private:
	T* a;
    };

    template <typename T>
    class arg<T const*>
    {
    public:
	arg(vector<T> const& a) : a(a.data()) {}
	arg(T const* a) : a(a) {}

	/**
	 * returns const pointer to CUDA device memory array
	 */
	T const* operator*()
	{
	    return a;
	}

    private:
	T const* a;
    };

    #endif /* ! __CUDACC__ */

    template <
	typename T,
	typename U = void(),
	typename V = void(),
	typename W = void()
	>
    class function;

    } // namespace cuda


    #define BOOST_PP_FILENAME_1 "cuda_wrapper/function.hpp"
    #define BOOST_PP_ITERATION_LIMITS (1, CUDA_FUNCTION_MAX_ARGS)
    #include BOOST_PP_ITERATE()

    #endif /* ! CUDA_FUNCTION_HPP */

#elif BOOST_PP_ITERATION_DEPTH() == 1

    #define CUDA_FUNCTION_ARGS_1 BOOST_PP_FRAME_ITERATION(1)

    namespace cuda
    {

    /**
     * CUDA kernel execution wrapper for n-ary device function
     */
    template <
	BOOST_PP_ENUM_PARAMS(CUDA_FUNCTION_ARGS_1, typename T)
	>
    class function<
	void (BOOST_PP_ENUM_PARAMS(CUDA_FUNCTION_ARGS_1, T)),
	void(),
	void(),
	void()
	>
    {
    public:
	typedef void T (BOOST_PP_ENUM_PARAMS(CUDA_FUNCTION_ARGS_1, T));

    public:
	function(T *f) : f(f) {}

    #ifndef __CUDACC__

	/**
	 * execute kernel
	 */
	void operator()(BOOST_PP_ENUM_BINARY_PARAMS(CUDA_FUNCTION_ARGS_1, arg<T, > x))
	{
	    // properly align CUDA device function arguments
	    struct {
		#define DECL_ARG(z, n, x) T##n x##n;
		BOOST_PP_REPEAT(CUDA_FUNCTION_ARGS_1, DECL_ARG, a)
		#undef DECL_ARG
	    } args = {
		BOOST_PP_ENUM_PARAMS(CUDA_FUNCTION_ARGS_1, *x)
	    };
	    // push aligned arguments onto CUDA execution stack
	    CUDA_CALL(cudaSetupArgument(&args, sizeof(args), 0));
	    // launch CUDA device function
	    CUDA_CALL(cudaLaunch(reinterpret_cast<const char *>(f)));
	}

    #endif /* ! __CUDACC__ */

    private:
	T *f;
    };


    } // namespace cuda

    #define BOOST_PP_FILENAME_2 "cuda_wrapper/function.hpp"
    #if CUDA_FUNCTION_VARY_ARGS
    #define BOOST_PP_ITERATION_LIMITS (1, CUDA_FUNCTION_MAX_ARGS)
    #else
    #define BOOST_PP_ITERATION_LIMITS (CUDA_FUNCTION_ARGS_1, CUDA_FUNCTION_ARGS_1)
    #endif
    #include BOOST_PP_ITERATE()

    #undef CUDA_FUNCTION_ARGS_1

#elif BOOST_PP_ITERATION_DEPTH() == 2

    #define CUDA_FUNCTION_ARGS_2 BOOST_PP_FRAME_ITERATION(2)

    namespace cuda
    {

    template <
	BOOST_PP_ENUM_PARAMS(CUDA_FUNCTION_ARGS_1, typename T),
	BOOST_PP_ENUM_PARAMS(CUDA_FUNCTION_ARGS_2, typename U)
	>
    class function<
	void (BOOST_PP_ENUM_PARAMS(CUDA_FUNCTION_ARGS_1, T)),
	void (BOOST_PP_ENUM_PARAMS(CUDA_FUNCTION_ARGS_2, U)),
	void (),
	void ()
	>
    {
    public:
	typedef void T (BOOST_PP_ENUM_PARAMS(CUDA_FUNCTION_ARGS_1, T));
	typedef void U (BOOST_PP_ENUM_PARAMS(CUDA_FUNCTION_ARGS_2, U));

    public:
	function(T *f1, U *f2) : f1(f1), f2(f2) {}

    #ifndef __CUDACC__
	void operator()(BOOST_PP_ENUM_BINARY_PARAMS(CUDA_FUNCTION_ARGS_1, arg<T, > x))
	{
	    f1(BOOST_PP_ENUM_PARAMS(CUDA_FUNCTION_ARGS_1, x));
	}
	void operator()(BOOST_PP_ENUM_BINARY_PARAMS(CUDA_FUNCTION_ARGS_2, arg<U, > x))
	{
	    f2(BOOST_PP_ENUM_PARAMS(CUDA_FUNCTION_ARGS_2, x));
	}
    #endif

    private:
	function<T> f1;
	function<U> f2;
    };

    } // namespace cuda

    #define BOOST_PP_FILENAME_3 "cuda_wrapper/function.hpp"
    #if CUDA_FUNCTION_VARY_ARGS
    #define BOOST_PP_ITERATION_LIMITS (1, CUDA_FUNCTION_MAX_ARGS)
    #else
    #define BOOST_PP_ITERATION_LIMITS (CUDA_FUNCTION_ARGS_2, CUDA_FUNCTION_ARGS_2)
    #endif
    #include BOOST_PP_ITERATE()

    #undef CUDA_FUNCTION_ARGS_2

#elif BOOST_PP_ITERATION_DEPTH() == 3

    #define CUDA_FUNCTION_ARGS_3 BOOST_PP_FRAME_ITERATION(3)

    namespace cuda
    {

    template <
	BOOST_PP_ENUM_PARAMS(CUDA_FUNCTION_ARGS_1, typename T),
	BOOST_PP_ENUM_PARAMS(CUDA_FUNCTION_ARGS_2, typename U),
	BOOST_PP_ENUM_PARAMS(CUDA_FUNCTION_ARGS_3, typename V)
	>
    class function<
	void (BOOST_PP_ENUM_PARAMS(CUDA_FUNCTION_ARGS_1, T)),
	void (BOOST_PP_ENUM_PARAMS(CUDA_FUNCTION_ARGS_2, U)),
	void (BOOST_PP_ENUM_PARAMS(CUDA_FUNCTION_ARGS_3, V)),
	void ()
	>
    {
    public:
	typedef void T (BOOST_PP_ENUM_PARAMS(CUDA_FUNCTION_ARGS_1, T));
	typedef void U (BOOST_PP_ENUM_PARAMS(CUDA_FUNCTION_ARGS_2, U));
	typedef void V (BOOST_PP_ENUM_PARAMS(CUDA_FUNCTION_ARGS_3, V));

    public:
	function(T *f1, U *f2, V *f3) : f1(f1), f2(f2), f3(f3) {}

    #ifndef __CUDACC__
	void operator()(BOOST_PP_ENUM_BINARY_PARAMS(CUDA_FUNCTION_ARGS_1, arg<T, > x))
	{
	    f1(BOOST_PP_ENUM_PARAMS(CUDA_FUNCTION_ARGS_1, x));
	}
	void operator()(BOOST_PP_ENUM_BINARY_PARAMS(CUDA_FUNCTION_ARGS_2, arg<U, > x))
	{
	    f2(BOOST_PP_ENUM_PARAMS(CUDA_FUNCTION_ARGS_2, x));
	}
	void operator()(BOOST_PP_ENUM_BINARY_PARAMS(CUDA_FUNCTION_ARGS_3, arg<V, > x))
	{
	    f3(BOOST_PP_ENUM_PARAMS(CUDA_FUNCTION_ARGS_3, x));
	}
    #endif

    private:
	function<T> f1;
	function<U> f2;
	function<V> f3;
    };

    } // namespace cuda

    #define BOOST_PP_FILENAME_4 "cuda_wrapper/function.hpp"
    #if CUDA_FUNCTION_VARY_ARGS
    #define BOOST_PP_ITERATION_LIMITS (1, CUDA_FUNCTION_MAX_ARGS)
    #else
    #define BOOST_PP_ITERATION_LIMITS (CUDA_FUNCTION_ARGS_3, CUDA_FUNCTION_ARGS_3)
    #endif
    #include BOOST_PP_ITERATE()

    #undef CUDA_FUNCTION_ARGS_3

#elif BOOST_PP_ITERATION_DEPTH() == 4

    #define CUDA_FUNCTION_ARGS_4 BOOST_PP_FRAME_ITERATION(4)

    namespace cuda
    {

    template <
	BOOST_PP_ENUM_PARAMS(CUDA_FUNCTION_ARGS_1, typename T),
	BOOST_PP_ENUM_PARAMS(CUDA_FUNCTION_ARGS_2, typename U),
	BOOST_PP_ENUM_PARAMS(CUDA_FUNCTION_ARGS_3, typename V),
	BOOST_PP_ENUM_PARAMS(CUDA_FUNCTION_ARGS_4, typename W)
	>
    class function<
	void (BOOST_PP_ENUM_PARAMS(CUDA_FUNCTION_ARGS_1, T)),
	void (BOOST_PP_ENUM_PARAMS(CUDA_FUNCTION_ARGS_2, U)),
	void (BOOST_PP_ENUM_PARAMS(CUDA_FUNCTION_ARGS_3, V)),
	void (BOOST_PP_ENUM_PARAMS(CUDA_FUNCTION_ARGS_4, W))
	>
    {
    public:
	typedef void T (BOOST_PP_ENUM_PARAMS(CUDA_FUNCTION_ARGS_1, T));
	typedef void U (BOOST_PP_ENUM_PARAMS(CUDA_FUNCTION_ARGS_2, U));
	typedef void V (BOOST_PP_ENUM_PARAMS(CUDA_FUNCTION_ARGS_3, V));
	typedef void W (BOOST_PP_ENUM_PARAMS(CUDA_FUNCTION_ARGS_4, W));

    public:
	function(T *f1, U *f2, V *f3, W *f4) : f1(f1), f2(f2), f3(f3), f4(f4) {}

    #ifndef __CUDACC__
	void operator()(BOOST_PP_ENUM_BINARY_PARAMS(CUDA_FUNCTION_ARGS_1, arg<T, > x))
	{
	    f1(BOOST_PP_ENUM_PARAMS(CUDA_FUNCTION_ARGS_1, x));
	}
	void operator()(BOOST_PP_ENUM_BINARY_PARAMS(CUDA_FUNCTION_ARGS_2, arg<U, > x))
	{
	    f2(BOOST_PP_ENUM_PARAMS(CUDA_FUNCTION_ARGS_2, x));
	}
	void operator()(BOOST_PP_ENUM_BINARY_PARAMS(CUDA_FUNCTION_ARGS_3, arg<V, > x))
	{
	    f3(BOOST_PP_ENUM_PARAMS(CUDA_FUNCTION_ARGS_3, x));
	}
	void operator()(BOOST_PP_ENUM_BINARY_PARAMS(CUDA_FUNCTION_ARGS_4, arg<W, > x))
	{
	    f4(BOOST_PP_ENUM_PARAMS(CUDA_FUNCTION_ARGS_4, x));
	}
    #endif

    private:
	function<T> f1;
	function<U> f2;
	function<V> f3;
	function<W> f4;
    };

    } // namespace cuda

    #undef CUDA_FUNCTION_ARGS_4

#endif /* !BOOST_PP_IS_ITERATING */
