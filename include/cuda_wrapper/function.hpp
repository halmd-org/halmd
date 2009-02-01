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
    #include <boost/tuple/tuple.hpp>
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

    #endif /* ! __CUDACC__ */

    template <typename T0,
	      typename T1 = void,
	      typename T2 = void,
	      typename T3 = void,
	      typename T4 = void>
    class function;

    template <typename T0, typename T1>
    class function<T0, T1> :
	function<T0>, function<T1>
    {
	typedef function<T0> F0;
	typedef function<T1> F1;

    public:
	function(T0* f0, T1* f1) :
	    F0(f0), F1(f1) {}
    #ifndef __CUDACC__
	using F0::operator();
	using F1::operator();
    #endif
    };

    template <typename T0, typename T1, typename T2>
    class function<T0, T1, T2> :
	function<T0>, function<T1>, function<T2>
    {
	typedef function<T0> F0;
	typedef function<T1> F1;
	typedef function<T2> F2;

    public:
	function(T0* f0, T1* f1, T2* f2) :
	    F0(f0), F1(f1), F2(f2) {}
    #ifndef __CUDACC__
	using F0::operator();
	using F1::operator();
	using F2::operator();
    #endif
    };

    template <typename T0, typename T1, typename T2, typename T3>
    class function<T0, T1, T2, T3> :
	function<T0>, function<T1>, function<T2>, function<T3>
    {
	typedef function<T0> F0;
	typedef function<T1> F1;
	typedef function<T2> F2;
	typedef function<T3> F3;

    public:
	function(T0* f0, T1* f1, T2* f2, T3* f3) :
	    F0(f0), F1(f1), F2(f2), F3(f3) {}
    #ifndef __CUDACC__
	using F0::operator();
	using F1::operator();
	using F2::operator();
	using F3::operator();
    #endif
    };


    template <typename T0, typename T1, typename T2, typename T3, typename T4>
    class function :
	function<T0>, function<T1>, function<T2>, function<T3>, function<T4>
    {
	typedef function<T0> F0;
	typedef function<T1> F1;
	typedef function<T2> F2;
	typedef function<T3> F3;
	typedef function<T4> F4;

    public:
	function(T0* f0, T1* f1, T2* f2, T3* f3, T4* f4) :
	    F0(f0), F1(f1), F2(f2), F3(f3), F4(f4) {}
    #ifndef __CUDACC__
	using F0::operator();
	using F1::operator();
	using F2::operator();
	using F3::operator();
	using F4::operator();
    #endif
    };

    } // namespace cuda

    #define BOOST_PP_FILENAME_1 <cuda_wrapper/function.hpp>
    #define BOOST_PP_ITERATION_LIMITS (1, CUDA_FUNCTION_MAX_ARGS)
    #include BOOST_PP_ITERATE()

    #endif /* ! CUDA_FUNCTION_HPP */

#elif BOOST_PP_ITERATION_DEPTH() == 1

    #define CUDA_FUNCTION_ARGS BOOST_PP_FRAME_ITERATION(1)
    #define CUDA_FUNCTION_PARAM_VALUES BOOST_PP_ENUM_PARAMS(CUDA_FUNCTION_ARGS, x)
    #define CUDA_FUNCTION_PARAM_TYPES BOOST_PP_ENUM_PARAMS(CUDA_FUNCTION_ARGS, T)

    namespace cuda
    {

    /**
     * CUDA kernel execution wrapper for n-ary device function
     */
    template <BOOST_PP_ENUM_PARAMS(CUDA_FUNCTION_ARGS, typename T)>
    class function<void (CUDA_FUNCTION_PARAM_TYPES)>
    {
    public:
	typedef void T (CUDA_FUNCTION_PARAM_TYPES);

    public:
	function(T *f) : f(f) {}

    #ifndef __CUDACC__

	/**
	 * execute kernel
	 */
	void operator()(BOOST_PP_ENUM_BINARY_PARAMS(CUDA_FUNCTION_ARGS, T, x))
	{
	    // properly align CUDA device function arguments
	    boost::tuple<CUDA_FUNCTION_PARAM_TYPES> args(CUDA_FUNCTION_PARAM_VALUES);
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

    #undef CUDA_FUNCTION_ARGS
    #undef CUDA_FUNCTION_PARAM_TYPES
    #undef CUDA_FUNCTION_PARAM_VALUES

#endif /* !BOOST_PP_IS_ITERATING */
