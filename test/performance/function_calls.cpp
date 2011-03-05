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

#define BOOST_TEST_MODULE function_calls
#include <boost/test/unit_test.hpp>

#include <boost/signals2.hpp>
#include <iomanip> // setw

#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>
#include <test/tools/lua.hpp>
#include <test/performance/function_calls_noop.hpp>
#ifdef WITH_CUDA
# include <test/performance/function_calls_noop_kernel.hpp>
#endif

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
     */
    printer(string const& method, size_t iterations)
      : iterations_(iterations)
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
        BOOST_CHECK( factor >= 1 );
    }

private:
    size_t iterations_;
};

// use macros to allow stringification for Lua commands
#define I1E7 10000000
// CUDA kernel launches are slow…
#define I1E6 1000000

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
 * Measure boost::function call
 */
BOOST_AUTO_TEST_CASE( boost_function )
{
    function<void (double)> f(bind_noop());
    printer p("boost::function", I1E7);
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
 * Measure vector<boost::function> call
 */
BOOST_AUTO_TEST_CASE( vector_boost_function )
{
    vector<function<void (double)> > f;
    bind_noop(f);
    printer p("vector<boost::function>", I1E7);
    // warm up
    for (size_t i = 0; i < I1E7; ++i) {
        vector<function<void (double)> >::const_iterator j, je = f.end();
        for (j = f.begin(); j != je; ++j) {
            (*j)(42.);
        }
    }
    scoped_timer<timer> timer(p);
    // benchmark
    for (size_t i = 0; i < I1E7; ++i) {
        vector<function<void (double)> >::const_iterator j, je = f.end();
        for (j = f.begin(); j != je; ++j) {
            (*j)(42.);
        }
    }
}

/**
 * Measure vector<boost::function> call with 10 slots
 */
BOOST_AUTO_TEST_CASE( vector_boost_function_10 )
{
    vector<function<void (double)> > f;
    for (size_t i = 0; i < 10; ++i) {
        bind_noop(f);
    }
    printer p("vector<boost::function> (10 slots)", I1E7);
    // warm up
    for (size_t i = 0; i < I1E7; ++i) {
        vector<function<void (double)> >::const_iterator j, je = f.end();
        for (j = f.begin(); j != je; ++j) {
            (*j)(42.);
        }
    }
    scoped_timer<timer> timer(p);
    // benchmark
    for (size_t i = 0; i < I1E7; ++i) {
        vector<function<void (double)> >::const_iterator j, je = f.end();
        for (j = f.begin(); j != je; ++j) {
            (*j)(42.);
        }
    }
}

/**
 * Measure boost::signals2 call
 */
BOOST_AUTO_TEST_CASE( boost_signals2 )
{
    using namespace boost::signals2;
    signal<void (double)> sig;
    sig.connect(bind_noop());
    printer p("boost::signals2", I1E7);
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
 * Measure boost::signals2 call with 10 slots
 */
BOOST_AUTO_TEST_CASE( boost_signals2_10 )
{
    using namespace boost::signals2;
    signal<void (double)> sig;
    for (size_t i = 0; i < 10; ++i) {
        sig.connect(bind_noop());
    }
    printer p("boost::signals2 (10 slots)", I1E6);
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

#ifdef WITH_CUDA

/**
 * Measure CUDA kernel call
 */
BOOST_AUTO_TEST_CASE( cuda_kernel_1_128 )
{
    printer p("CUDA kernel (1×128)", I1E6);
    // warm up
    for (size_t i = 0; i < I1E6; ++i) {
        launch_noop_kernel(1, 128, 42.);
    }
    scoped_timer<timer> timer(p);
    // benchmark
    for (size_t i = 0; i < I1E6; ++i) {
        launch_noop_kernel(1, 128, 42.);
    }
}

/**
 * Measure CUDA kernel call
 */
BOOST_AUTO_TEST_CASE( cuda_kernel_16_512 )
{
    printer p("CUDA kernel (16×512)", I1E6);
    // warm up
    for (size_t i = 0; i < I1E6; ++i) {
        launch_noop_kernel(16, 512, 42.);
    }
    scoped_timer<timer> timer(p);
    // benchmark
    for (size_t i = 0; i < I1E6; ++i) {
        launch_noop_kernel(16, 512, 42.);
    }
}

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
    scoped_timer<timer> timer(p);
    // benchmark
    for (size_t i = 0; i < I1E6; ++i) {
        cuda::configure(1, 128);
        noop_kernel(42.);
    }
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
    scoped_timer<timer> timer(p);
    // benchmark
    for (size_t i = 0; i < I1E6; ++i) {
        cuda::configure(16, 512);
        noop_kernel(42.);
    }
}

#endif /* WITH_CUDA */

/**
 * Measure global Lua function call
 */
BOOST_FIXTURE_TEST_CASE( global_lua_function, lua_test_fixture )
{
    LUA_CHECK( "noop = function() end" );
    printer p("Lua function (global)", I1E7);
    // warm up
    LUA_CHECK( "local noop = noop for i = 1, " xstr(I1E7) " do noop(42.) end" );
    scoped_timer<timer> timer(p);
    // benchmark
    LUA_CHECK( "for i = 1, " xstr(I1E7) " do noop(42.) end" );
}

/**
 * Measure local Lua function call
 */
BOOST_FIXTURE_TEST_CASE( local_lua_function, lua_test_fixture )
{
    LUA_CHECK( "noop = function() end" );
    printer p("Lua function (local)", I1E7);
    // warm up
    LUA_CHECK( "local noop = noop for i = 1, " xstr(I1E7) " do noop(42.) end" );
    scoped_timer<timer> timer(p);
    // benchmark
    LUA_CHECK( "local noop = noop for i = 1, " xstr(I1E7) " do noop(42.) end" );
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
 * Measure Luabind function call
 */
BOOST_FIXTURE_TEST_CASE( luabind_function, lua_test_fixture )
{
    using namespace luabind;
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
