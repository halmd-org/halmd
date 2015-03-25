/*
 * Copyright © 2011  Peter Colberg
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

#include <halmd/config.hpp>

#define BOOST_TEST_MODULE function_calls
#include <boost/test/unit_test.hpp>

#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/signal.hpp>
#include <halmd/utility/timer.hpp>
#include <test/tools/ctest.hpp>
#include <test/tools/lua.hpp>
#include <test/performance/function_calls_extern.hpp>
#ifdef HALMD_WITH_GPU
# include <test/tools/cuda.hpp>
# include <test/performance/function_calls_extern_kernel.hpp>
#endif

#include <iomanip>

/**
 * This test suite compares the invokation times of different methods of
 * calling an empty function: Directly, via a functor or via Luabind, or
 * as a CUDA kernel. The goal is to give an idea of the overhead of the
 * different methods over a direct C++ function call.
 */

// http://gcc.gnu.org/onlinedocs/cpp/Stringification.html
#define xstr(s) str(s)
#define str(s) #s

using namespace boost;
using namespace halmd; // scoped_timer, timer
using namespace std;

/**
 * Print result, and compare to direct function call.
 */
class printer
{
public:
    /**
     * Construct result printer for use with scoped_timer.
     *
     * @param method name of function call method
     * @param iterations number of function calls
     * @param check perform test if the performance test was slower than the C++ test
     */
    printer(string const& method, size_t iterations, bool check = true)
      : iterations_(iterations)
      , check_(check)
    {
        BOOST_TEST_MESSAGE( method );
    }

    /**
     * Invoked when scoped_timer exits scope.
     *
     * @param elapsed elapsed time in seconds
     */
    void operator()(double elapsed) const
    {
        // time in seconds per function call
        double result = elapsed / iterations_;
        // store result time of C++ function call
        static double first = result;
        // result time relative to C++ function call
        double factor = result / first;

        BOOST_TEST_MESSAGE( setprecision(2) << fixed << setw(15) << result * 1.e9 << " ns" << setw(15) << factor << "" );

        // C++ function call should be fastest
        if (check_) {
            BOOST_CHECK( factor >= 1 );
        }
    }

private:
    size_t iterations_;
    bool check_;
};

// use macros to allow stringification for Lua commands
#define I1E7 10000000
// CUDA kernel launches are slow…
#define I1E6 1000000
// synchroneous CUDA kernel launches are even slower…
#define I1E5 100000

BOOST_AUTO_TEST_SUITE( host )

/**
 * Measure direct C++ function call
 */
BOOST_AUTO_TEST_CASE( cpp_function )
{
    printer p("C++ function", I1E7);
    // warm up
    for (size_t i = 0; i < I1E7; ++i) {
        noop(42.);
    }
    scoped_timer<timer> timer(p);
    // benchmark
    for (size_t i = 0; i < I1E7; ++i) {
        noop(42.);
    }
}

/**
 * Measure std::function call
 */
BOOST_AUTO_TEST_CASE( std_function )
{
    std::function<void (double)> f(bind_noop());
    printer p("std::function", I1E7);
    // warm up
    for (size_t i = 0; i < I1E7; ++i) {
        f(42.);
    }
    scoped_timer<timer> timer(p);
    // benchmark
    for (size_t i = 0; i < I1E7; ++i) {
        f(42.);
    }
}

/**
 * Measure halmd::signal call
 */
BOOST_AUTO_TEST_CASE( halmd_signal )
{
    halmd::signal<void (double)> sig;
    bind_noop(sig);
    printer p("halmd::signal", I1E7);
    // warm up
    for (size_t i = 0; i < I1E7; ++i) {
        sig(42.);
    }
    scoped_timer<timer> timer(p);
    // benchmark
    for (size_t i = 0; i < I1E7; ++i) {
        sig(42.);
    }
}

/**
 * Measure halmd::signal call with 10 slots
 */
BOOST_AUTO_TEST_CASE( halmd_signal_10 )
{
    halmd::signal<void (double)> sig;
    for (size_t i = 0; i < 10; ++i) {
        bind_noop(sig);
    }
    printer p("halmd::signal (10 slots)", I1E6);
    // warm up
    for (size_t i = 0; i < I1E6; ++i) {
        sig(42.);
    }
    scoped_timer<timer> timer(p);
    // benchmark
    for (size_t i = 0; i < I1E6; ++i) {
        sig(42.);
    }
}

// see Ch. 26.1 in Programming in Lua by R. Ierusalimschy
// http://www.lua.org/pil/26.1.html
static int lua_noop(lua_State* L)
{
    double dummy = lua_tonumber(L, 1);  // get argument
    noop(dummy);
    return 0;                           // number of results
}

/**
 * Measure Lua C function call
 */
BOOST_FIXTURE_TEST_CASE( lua_cfunction, lua_test_fixture )
{
    lua_pushcfunction(L, lua_noop);
    lua_setglobal(L, "noop");
    printer p("Lua C function", I1E7);
    // warm up
    LUA_CHECK( "local noop = noop for i = 1, " xstr(I1E7) " do noop(42.) end" );
    scoped_timer<timer> timer(p);
    // benchmark
    LUA_CHECK( "local noop = noop for i = 1, " xstr(I1E7) " do noop(42.) end" );
}

/**
 * Measure local Lua function call
 */
BOOST_FIXTURE_TEST_CASE( local_lua_function, lua_test_fixture )
{
    LUA_CHECK( "noop = function() end" );
    printer p("Lua function (local)", I1E7, false); // LuaJIT may generate faster code than C++
    // warm up
    LUA_CHECK( "local noop = noop for i = 1, " xstr(I1E7) " do noop(42.) end" );
    scoped_timer<timer> timer(p);
    // benchmark
    LUA_CHECK( "local noop = noop for i = 1, " xstr(I1E7) " do noop(42.) end" );
}

/**
 * Measure global Lua function call
 */
BOOST_FIXTURE_TEST_CASE( global_lua_function, lua_test_fixture )
{
    LUA_CHECK( "noop = function() end" );
    printer p("Lua function (global)", I1E7, false); // LuaJIT may generate faster code than C++
    // warm up
    LUA_CHECK( "for i = 1, " xstr(I1E7) " do noop(42.) end" );
    scoped_timer<timer> timer(p);
    // benchmark
    LUA_CHECK( "for i = 1, " xstr(I1E7) " do noop(42.) end" );
}

/**
 * Measure Luabind function call
 */
BOOST_FIXTURE_TEST_CASE( luaponte_function, lua_test_fixture )
{
    using namespace luaponte;
    module(L)
    [
        def("noop", &noop)
    ];
    printer p("Luabind function", I1E7);
    // warm up
    LUA_CHECK( "local noop = noop for i = 1, " xstr(I1E7) " do noop(42.) end" );
    scoped_timer<timer> timer(p);
    // benchmark
    LUA_CHECK( "local noop = noop for i = 1, " xstr(I1E7) " do noop(42.) end" );
}

BOOST_AUTO_TEST_SUITE_END() // host

#ifdef HALMD_WITH_GPU

BOOST_FIXTURE_TEST_SUITE( gpu, set_cuda_device )

/**
 * Measure cuda::function call
 */
BOOST_AUTO_TEST_CASE( cuda_function_1_128 )
{
    printer p("cuda::function (1×128)", I1E6);
    // warm up
    for (size_t i = 0; i < I1E6; ++i) {
        cuda::configure(1, 128);
        noop_kernel(42.);
    }
    cuda::thread::synchronize();
    scoped_timer<timer> timer(p);
    // benchmark
    for (size_t i = 0; i < I1E6; ++i) {
        cuda::configure(1, 128);
        noop_kernel(42.);
    }
    cuda::thread::synchronize();
}

/**
 * Measure cuda::function call
 */
BOOST_AUTO_TEST_CASE( cuda_function_16_512 )
{
    printer p("cuda::function (16×512)", I1E6);
    // warm up
    for (size_t i = 0; i < I1E6; ++i) {
        cuda::configure(16, 512);
        noop_kernel(42.);
    }
    cuda::thread::synchronize();
    scoped_timer<timer> timer(p);
    // benchmark
    for (size_t i = 0; i < I1E6; ++i) {
        cuda::configure(16, 512);
        noop_kernel(42.);
    }
    cuda::thread::synchronize();
}

/**
 * Measure cuda::function call
 */
BOOST_AUTO_TEST_CASE( cuda_function_16_512_synchronize )
{
    printer p("cuda::function (16×512) + cuda::thread::synchronize", I1E5);
    // warm up
    for (size_t i = 0; i < I1E5; ++i) {
        cuda::configure(16, 512);
        noop_kernel(42.);
        cuda::thread::synchronize();
    }
    scoped_timer<timer> timer(p);
    // benchmark
    for (size_t i = 0; i < I1E5; ++i) {
        cuda::configure(16, 512);
        noop_kernel(42.);
        cuda::thread::synchronize();
    }
}

/**
 * Measure CUDA kernel call
 */
BOOST_AUTO_TEST_CASE( cuda_kernel_1_128 )
{
    printer p("CUDA <<< >>> (1×128)", I1E6);
    // warm up
    for (size_t i = 0; i < I1E6; ++i) {
        launch_noop_kernel(1, 128, 42.);
    }
    cudaThreadSynchronize();
    scoped_timer<timer> timer(p);
    // benchmark
    for (size_t i = 0; i < I1E6; ++i) {
        launch_noop_kernel(1, 128, 42.);
    }
    cudaThreadSynchronize();
}

/**
 * Measure CUDA kernel call
 */
BOOST_AUTO_TEST_CASE( cuda_kernel_16_512 )
{
    printer p("CUDA <<< >>> (16×512)", I1E6);
    // warm up
    for (size_t i = 0; i < I1E6; ++i) {
        launch_noop_kernel(16, 512, 42.);
    }
    cudaThreadSynchronize();
    scoped_timer<timer> timer(p);
    // benchmark
    for (size_t i = 0; i < I1E6; ++i) {
        launch_noop_kernel(16, 512, 42.);
    }
    cudaThreadSynchronize();
}

/**
 * Measure CUDA kernel call
 */
BOOST_AUTO_TEST_CASE( cuda_kernel_16_512_synchronize )
{
    printer p("CUDA <<< >>> (16×512) + cudaThreadSynchronize", I1E5);
    // warm up
    for (size_t i = 0; i < I1E5; ++i) {
        launch_noop_kernel(16, 512, 42.);
        cudaThreadSynchronize();
    }
    scoped_timer<timer> timer(p);
    // benchmark
    for (size_t i = 0; i < I1E5; ++i) {
        launch_noop_kernel(16, 512, 42.);
        cudaThreadSynchronize();
    }
}

BOOST_AUTO_TEST_SUITE_END() // gpu

#endif /* HALMD_WITH_GPU */
