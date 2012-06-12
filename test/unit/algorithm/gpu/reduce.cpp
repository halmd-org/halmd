/*
 * Copyright © 2010-2012 Peter Colberg
 * Copyright © 2011 Felix Höfling
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

#define BOOST_TEST_MODULE reduce
#include <boost/test/unit_test.hpp>

#include <halmd/algorithm/gpu/reduce.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>
#include <halmd/utility/timer.hpp>
#include <test/tools/ctest.hpp>
#include <test/tools/cuda.hpp>
#include <test/tools/init.hpp>
#include <test/unit/algorithm/gpu/reduce_kernel.hpp>

#include <boost/iterator/counting_iterator.hpp>
#include <boost/test/parameterized_test.hpp>

#include <algorithm>
#include <array>
#include <functional>

BOOST_GLOBAL_FIXTURE( set_cuda_device )

/**
 * Compute sum of natural numbers using a unary reduction.
 */
template <typename accumulator_type>
static void reduce_sum()
{
    cuda::host::vector<float> h_v(
        boost::make_counting_iterator(1)
      , boost::make_counting_iterator(12345679)
    );
    BOOST_TEST_MESSAGE( "  summation of " << h_v.size() << " floats" );
    cuda::vector<float> g_v(h_v.size());
    BOOST_CHECK( cuda::copy(h_v.begin(), h_v.end(), g_v.begin()) == g_v.end());
    accumulator_type acc = halmd::reduce(
        &*g_v.begin()
      , &*g_v.end()
      , accumulator_type(0)
    );
    BOOST_CHECK_EQUAL(std::int64_t(double(acc())), 12345678LL * 12345679 / 2);
}

/**
 * Test »32 bit integer arithmetic using double-single floating point (48 bit mantissa).
 */
BOOST_AUTO_TEST_CASE( reduce_sum_to_double_single )
{
    reduce_sum<sum<float, halmd::dsfloat> >();
}

#ifdef HALMD_GPU_DOUBLE_PRECISION
/**
 * Test »32 bit integer arithmetic using double-precision floating point (53 bit mantissa).
 */
BOOST_AUTO_TEST_CASE( reduce_sum_to_double )
{
    reduce_sum<sum<float, double> >();
}
#endif

/**
 * Compute sum of natural numbers using a unary reduction.
 */
template <typename accumulator_type>
static void reduce_sum_with_constant()
{
    cuda::host::vector<float> h_v(
        boost::make_counting_iterator(1)
      , boost::make_counting_iterator(12345679)
    );
    BOOST_TEST_MESSAGE( "  summation of " << h_v.size() << " floats" );
    cuda::vector<float> g_v(h_v.size());
    BOOST_CHECK( cuda::copy(h_v.begin(), h_v.end(), g_v.begin()) == g_v.end());
    accumulator_type acc = halmd::reduce(
        std::make_tuple(&*g_v.begin(), -1)
      , std::make_tuple(&*g_v.end())
      , accumulator_type(0)
    );
    BOOST_CHECK_EQUAL(std::int64_t(double(acc())), -12345678LL * 12345679 / 2);
}

/**
 * Test »32 bit integer arithmetic using double-single floating point (48 bit mantissa).
 */
BOOST_AUTO_TEST_CASE( reduce_sum_with_constant_to_double_single )
{
    reduce_sum_with_constant<sum_with_constant<float, halmd::dsfloat> >();
}

#ifdef HALMD_GPU_DOUBLE_PRECISION
/**
 * Test »32 bit integer arithmetic using double-precision floating point (53 bit mantissa).
 */
BOOST_AUTO_TEST_CASE( reduce_sum_with_constant_to_double )
{
    reduce_sum_with_constant<sum_with_constant<float, double> >();
}
#endif

/**
 * Compute sum of squares of natural numbers using a binary reduction.
 */
template <typename accumulator_type>
static void reduce_sum_of_squares()
{
    cuda::host::vector<float> h_v1(
        boost::make_counting_iterator(1)
      , boost::make_counting_iterator(123457)
    );
    BOOST_TEST_MESSAGE( "  summation of squares of " << h_v1.size() << " floats" );
    cuda::host::vector<float> h_v2(h_v1.size());
    std::transform(
        h_v1.begin()
      , h_v1.end()
      , h_v2.begin()
      , [](float value) {
            return -value;
        }
    );
    cuda::vector<float> g_v1(h_v1.size());
    cuda::vector<float> g_v2(h_v2.size());
    BOOST_CHECK( cuda::copy(h_v1.begin(), h_v1.end(), g_v1.begin()) == g_v1.end());
    BOOST_CHECK( cuda::copy(h_v2.begin(), h_v2.end(), g_v2.begin()) == g_v2.end());
    accumulator_type acc = halmd::reduce(
        std::make_tuple(&*g_v1.begin(), &*g_v2.begin())
      , std::make_tuple(&*g_v1.end())
      , accumulator_type(0)
    );
    BOOST_CHECK_EQUAL(std::int64_t(double(acc())), -123456LL * 123457 * (2 * 123456 + 1) / 6);
}

/**
 * Test »32 bit integer arithmetic using double-single floating point (48 bit mantissa).
 */
BOOST_AUTO_TEST_CASE( reduce_sum_of_squares_to_double_single )
{
    reduce_sum_of_squares<sum_of_squares<float, halmd::dsfloat> >();
}

#ifdef HALMD_GPU_DOUBLE_PRECISION
/**
 * Test »32 bit integer arithmetic using double-precision floating point (53 bit mantissa).
 */
BOOST_AUTO_TEST_CASE( reduce_sum_of_squares_to_double )
{
    reduce_sum_of_squares<sum_of_squares<float, double> >();
}
#endif

/**
 * Compute sum of squares of natural numbers using a binary reduction.
 */
template <typename accumulator_type>
static void reduce_sum_of_squares_with_constant()
{
    cuda::host::vector<float> h_v1(
        boost::make_counting_iterator(1)
      , boost::make_counting_iterator(123457)
    );
    BOOST_TEST_MESSAGE( "  summation of squares of " << h_v1.size() << " floats" );
    cuda::host::vector<float> h_v2(h_v1.size());
    std::transform(
        h_v1.begin()
      , h_v1.end()
      , h_v2.begin()
      , [](float value) {
            return -value;
        }
    );
    cuda::vector<float> g_v1(h_v1.size());
    cuda::vector<float> g_v2(h_v2.size());
    BOOST_CHECK( cuda::copy(h_v1.begin(), h_v1.end(), g_v1.begin()) == g_v1.end());
    BOOST_CHECK( cuda::copy(h_v2.begin(), h_v2.end(), g_v2.begin()) == g_v2.end());
    accumulator_type acc = halmd::reduce(
        std::make_tuple(&*g_v1.begin(), &*g_v2.begin(), 7)
      , std::make_tuple(&*g_v1.end())
      , accumulator_type(0)
    );
    BOOST_CHECK_EQUAL(std::int64_t(double(acc())), -7 * 123456LL * 123457 * (2 * 123456 + 1) / 6);
}

/**
 * Test »32 bit integer arithmetic using double-single floating point (48 bit mantissa).
 */
BOOST_AUTO_TEST_CASE( reduce_sum_of_squares_with_constant_to_double_single )
{
    reduce_sum_of_squares_with_constant<sum_of_squares_with_constant<float, halmd::dsfloat> >();
}

#ifdef HALMD_GPU_DOUBLE_PRECISION
/**
 * Test »32 bit integer arithmetic using double-precision floating point (53 bit mantissa).
 */
BOOST_AUTO_TEST_CASE( reduce_sum_of_squares_with_constant_to_double )
{
    reduce_sum_of_squares_with_constant<sum_of_squares_with_constant<float, double> >();
}
#endif

/**
 * Compute sum of cubes of natural numbers using a binary reduction.
 */
template <typename accumulator_type>
static void reduce_sum_of_cubes()
{
    cuda::host::vector<float> h_v1(
        boost::make_counting_iterator(1)
      , boost::make_counting_iterator(1235)
    );
    BOOST_TEST_MESSAGE( "  summation of cubes of " << h_v1.size() << " floats" );
    cuda::host::vector<float> h_v2(h_v1.size());
    std::transform(
        h_v1.begin()
      , h_v1.end()
      , h_v2.begin()
      , [](float value) {
            return -value;
        }
    );
    cuda::host::vector<float> h_v3(h_v1.size());
    std::transform(
        h_v1.begin()
      , h_v1.end()
      , h_v3.begin()
      , [](float value) {
            return 7 * value;
        }
    );
    cuda::vector<float> g_v1(h_v1.size());
    cuda::vector<float> g_v2(h_v2.size());
    cuda::vector<float> g_v3(h_v3.size());
    BOOST_CHECK( cuda::copy(h_v1.begin(), h_v1.end(), g_v1.begin()) == g_v1.end());
    BOOST_CHECK( cuda::copy(h_v2.begin(), h_v2.end(), g_v2.begin()) == g_v2.end());
    BOOST_CHECK( cuda::copy(h_v3.begin(), h_v3.end(), g_v3.begin()) == g_v3.end());
    accumulator_type acc = halmd::reduce(
        std::make_tuple(&*g_v1.begin(), &*g_v2.begin(), &*g_v3.begin())
      , std::make_tuple(&*g_v1.end())
      , accumulator_type(0)
    );
    BOOST_CHECK_EQUAL(std::int64_t(double(acc())), -7 * 1234LL * 1234LL * 1235LL * 1235LL / 4);
}

/**
 * Test »32 bit integer arithmetic using double-single floating point (48 bit mantissa).
 */
BOOST_AUTO_TEST_CASE( reduce_sum_of_cubes_to_double_single )
{
    reduce_sum_of_cubes<sum_of_cubes<float, halmd::dsfloat> >();
}

#ifdef HALMD_GPU_DOUBLE_PRECISION
/**
 * Test »32 bit integer arithmetic using double-precision floating point (53 bit mantissa).
 */
BOOST_AUTO_TEST_CASE( reduce_sum_of_cubes_to_double )
{
    reduce_sum_of_cubes<sum_of_cubes<float, double> >();
}
#endif

/**
 * benchmark reduce kernel
 */
static void performance(std::size_t size)
{
#ifdef HALMD_GPU_DOUBLE_PRECISION
    typedef sum<float, double> accumulator_type;
#else
    typedef sum<float, halmd::dsfloat> accumulator_type;
#endif

    const int count_local = 20;
    const int count_global = 100;
    const int blocks = 30;
    const int threads = 64 << 3;

    try {
        cuda::host::vector<float> h_v(boost::make_counting_iterator(std::size_t(1)), boost::make_counting_iterator(size));
        cuda::vector<float> g_v(h_v.size());
        BOOST_CHECK( cuda::copy(h_v.begin(), h_v.end(), g_v.begin()) == g_v.end());

        std::array<halmd::accumulator<double>, 2> elapsed;
        for (int i = 0; i < count_local; ++i) {
            halmd::timer timer;
            // initialise kernel and allocate internal memory
            halmd::reduction<accumulator_type> reduce(blocks, threads);
            accumulator_type sum_local = reduce(&*g_v.begin(), &*g_v.end());
            elapsed[0](timer.elapsed());
            BOOST_CHECK_EQUAL(std::uint64_t(double(sum_local())), std::uint64_t(size - 1) * size / 2);
        }

        // pre-initialise kernel
        halmd::reduction<accumulator_type> reduce(blocks, threads);
        for (int i = 0; i < count_global; ++i) {
            halmd::timer timer;
            accumulator_type sum_global = reduce(&*g_v.begin(), &*g_v.end());
            elapsed[1](timer.elapsed());
            BOOST_CHECK_EQUAL(std::uint64_t(double(sum_global())), std::uint64_t(size - 1) * size / 2);
        }

        BOOST_TEST_MESSAGE("  summation of " << size << " floats: "
            << mean(elapsed[0]) * 1e3 << " ± " << error_of_mean(elapsed[0]) * 1e3 << " ms (local), "
            << mean(elapsed[1]) * 1e3 << " ± " << error_of_mean(elapsed[1]) * 1e3 << " ms (global)"
        );
    }
    catch (cuda::error const& e) {
        BOOST_CHECK(e.err == cudaErrorMemoryAllocation);
        BOOST_TEST_MESSAGE("Unsufficient device memory. Skip test for array of " << size << " bytes");
    }
    catch (std::bad_alloc const& e) {
        BOOST_TEST_MESSAGE("Unsufficient host memory. Skip test for array of " << size << " bytes");
    }
}

HALMD_TEST_INIT( init_unit_test_suite )
{
    using namespace boost::unit_test::framework;

    std::vector<std::size_t> sizes;
    for (std::size_t s = 1; s <= (1L << 29); s <<= 1) {
        sizes.push_back(s);
    }

    master_test_suite().add(BOOST_PARAM_TEST_CASE(&performance, sizes.begin(), sizes.end()));
}
