/*
 * Copyright © 2010-2012  Felix Höfling and Peter Colberg
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
#include <boost/test/parameterized_test.hpp>

#include <algorithm> // std::transform
#include <boost/array.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <functional> // std::negate

#include <halmd/algorithm/gpu/reduce.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>
#include <halmd/utility/timer.hpp>
#include <test/tools/cuda.hpp>
#include <test/tools/init.hpp>
#include <test/unit/algorithm/gpu/reduce_kernel.hpp>

using namespace boost;
using namespace halmd;

BOOST_GLOBAL_FIXTURE( set_cuda_device );

/**
 * Compute sum of natural numbers using a unary reduction.
 */
template <typename accumulator_type>
static void unary_reduce_float()
{
    cuda::host::vector<float> h_v(make_counting_iterator(1), make_counting_iterator(12345679));
    cuda::vector<float> g_v(h_v.size());
    cuda::copy(h_v, g_v);
    accumulator_type acc = reduce(g_v, accumulator_type());
    BOOST_CHECK_EQUAL(int64_t(double(acc())), 12345678LL * 12345679 / 2);
}

/**
 * Test »32 bit integer arithmetic using double-single floating point (48 bit mantissa).
 */
BOOST_AUTO_TEST_CASE( unary_reduce_float_to_double_single )
{
    unary_reduce_float<sum<float, dsfloat> >();
}

#ifdef HALMD_GPU_DOUBLE_PRECISION
/**
 * Test »32 bit integer arithmetic using double-precision floating point (53 bit mantissa).
 */
BOOST_AUTO_TEST_CASE( unary_reduce_float_to_double )
{
    unary_reduce_float<sum<float, double> >();
}
#endif

/**
 * Compute sum of squares of natural numbers using a binary reduction.
 */
template <typename accumulator_type>
static void binary_reduce_float()
{
    cuda::host::vector<float> h_v1(make_counting_iterator(1), make_counting_iterator(123457));
    cuda::host::vector<float> h_v2(h_v1.size());
    std::transform(h_v1.begin(), h_v1.end(), h_v2.begin(), std::negate<float>());
    cuda::vector<float> g_v1(h_v1.size());
    cuda::vector<float> g_v2(h_v2.size());
    cuda::copy(h_v1, g_v1);
    cuda::copy(h_v2, g_v2);
    accumulator_type acc = reduce(g_v1, g_v2, accumulator_type());
    BOOST_CHECK_EQUAL(int64_t(double(acc())), -123456LL * 123457 * (2 * 123456 + 1) / 6);
}

/**
 * Test »32 bit integer arithmetic using double-single floating point (48 bit mantissa).
 */
BOOST_AUTO_TEST_CASE( binary_reduce_float_to_double_single )
{
    binary_reduce_float<sum_of_squares<float, dsfloat> >();
}

#ifdef HALMD_GPU_DOUBLE_PRECISION
/**
 * Test »32 bit integer arithmetic using double-precision floating point (53 bit mantissa).
 */
BOOST_AUTO_TEST_CASE( binary_reduce_float_to_double )
{
    binary_reduce_float<sum_of_squares<float, double> >();
}
#endif

/**
 * benchmark reduce kernel
 */
static void performance(size_t size)
{
#ifdef HALMD_GPU_DOUBLE_PRECISION
    typedef sum<float, double> accumulator_type;
#else
    typedef sum<float, dsfloat> accumulator_type;
#endif

    const int count_local = 20;
    const int count_global = 100;
    const int blocks = 30;
    const int threads = 64 << 3;

    try {
        cuda::host::vector<float> h_v(make_counting_iterator(size_t(1)), make_counting_iterator(size));
        cuda::vector<float> g_v(h_v.size());
        cuda::copy(h_v, g_v);

        array<accumulator<double>, 2> elapsed;
        for (int i = 0; i < count_local; ++i) {
            halmd::timer timer;
            // initialise kernel and allocate internal memory
            reduction<accumulator_type> reduce(blocks, threads);
            accumulator_type sum_local = reduce(g_v, accumulator_type());
            elapsed[0](timer.elapsed());
            BOOST_CHECK_EQUAL(uint64_t(double(sum_local())), uint64_t(size - 1) * size / 2);
        }

        // pre-initialise kernel
        reduction<accumulator_type> reduce(blocks, threads);
        for (int i = 0; i < count_global; ++i) {
            halmd::timer timer;
            accumulator_type sum_global = reduce(g_v, accumulator_type());
            elapsed[1](timer.elapsed());
            BOOST_CHECK_EQUAL(uint64_t(double(sum_global())), uint64_t(size - 1) * size / 2);
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

    std::vector<size_t> sizes;
    for (size_t s = 1; s <= (1L << 29); s <<= 1) {
        sizes.push_back(s);
    }

    master_test_suite().add(BOOST_PARAM_TEST_CASE(&performance, sizes.begin(), sizes.end()));
}
