/* test random number generators in combination with accumulator class
 *
 * Copyright (C) 2010  Felix Höfling
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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_random_distributions
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>

#include <cuda_wrapper.hpp>
#include <time.h>
#include <stdexcept>

#include <halmd/math/accum.hpp>
#include <halmd/rng/gsl_rng.hpp>
#include <halmd/rng/rand48.hpp>

const unsigned BLOCKS = 64;
const unsigned THREADS = 128;

// "global fixture:" select CUDA device
struct set_cuda_device {
    set_cuda_device() {
        try {
            cuda::device::set(0);
        }
        catch (cuda::error const& e) {
            std::cerr << "CUDA error: " << e.what() << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    ~set_cuda_device() {}  // release device here?
};

BOOST_GLOBAL_FIXTURE( set_cuda_device );

void test_accum()
{
    halmd::accumulator<double> a;

    for (unsigned i=0; i <= 1; i++) {
        for (unsigned j=0; j < 10; j++) {
            a += j;
        }

        BOOST_CHECK_EQUAL(a.count(), 10u);
        BOOST_CHECK_CLOSE_FRACTION(a.mean(), 4.5, 1e-14);
        BOOST_CHECK_CLOSE_FRACTION(a.var(), 8.25, 1e-14);
        BOOST_CHECK_CLOSE_FRACTION(a.std(), 2.8722813232690143, 1e-14);
        BOOST_CHECK_CLOSE_FRACTION(a.err(), 0.9574271077563381, 1e-14);

        a.clear();  // check clear() function
    }
}

void test_rand48_gpu( unsigned long count )
{
    unsigned seed = time(NULL);

    try {
        // seed GPU random number generator
        halmd::rand48 rng(cuda::config(BLOCKS, THREADS));
        rng.set(seed);
        cuda::thread::synchronize();

        // parallel GPU rand48
        cuda::vector<float> g_array(count);
        rng.uniform(g_array);
        cuda::thread::synchronize();

        cuda::host::vector<float> h_array(count);
        cuda::copy(g_array, h_array);

        halmd::accumulator<float> a;
        for (unsigned long i=0; i < count; i++) {
            a += h_array[i];
        }

        // mean = 1/2, std = 1/12
        // use tolerance = six sigma, so the test passes with 97% probability
        double tol = 6 * sqrt(1. / (count - 1) / (1./12));
        BOOST_CHECK_EQUAL(a.count(), count);
        BOOST_CHECK_CLOSE_FRACTION(a.mean(), 0.5, tol / .5);
        BOOST_CHECK_CLOSE_FRACTION(a.var(), 1./12, tol / (1./12));
    }
    catch (cuda::error const& e) {
        BOOST_FAIL("(CUDA error) " << e.what());
    }
    catch (std::exception const& e) {
        BOOST_FAIL(e.what());
    }
}

void test_gsl_rng( unsigned long count )
{
    unsigned seed = time(NULL);

    // seed GPU random number generator
    halmd::gsl::gfsr4 rng;
    rng.set(seed);

    // test uniform distribution
    halmd::accumulator<double> a;
    for (unsigned i=0; i < count; i++)
        a += rng.uniform();

    // mean = 1/2, std = 1/12
    // use tolerance = six sigma, so the test passes with 97% probability
    double tol = 6 * sqrt(1. / (count - 1) / (1./12));
    BOOST_CHECK_EQUAL(a.count(), count);
    BOOST_CHECK_CLOSE_FRACTION(a.mean(), .5, tol / .5);
    BOOST_CHECK_CLOSE_FRACTION(a.var(), 1./12, tol / (1./12));

    // test Gaussian distribution
    BOOST_REQUIRE_MESSAGE(count % 2 == 0, "require even number of samples");
    a.clear();
    for (unsigned i=0; i < count; i+=2) {
        double x1, x2;
        rng.gaussian(x1, x2, 1.);
        a += 1 + x1;
        a += 1 + x2;
    }

    // mean = 1, std = 1
    tol = 6 * sqrt(1. / (count - 1));  // tolerance = six sigma
    BOOST_CHECK_EQUAL(a.count(), count);
    BOOST_CHECK_CLOSE_FRACTION(a.mean(), 1, tol);
    BOOST_CHECK_CLOSE_FRACTION(a.var(), 1, tol);
}

int init_unit_test_suite()
{
    using namespace boost::unit_test::framework;

    std::vector<unsigned long> counts;
    counts.push_back(100);
    counts.push_back(1000);
    counts.push_back(10000);
    counts.push_back(100000);
    counts.push_back(1000000);
    counts.push_back(10000000);
    counts.push_back(100000000);

    master_test_suite().add(BOOST_TEST_CASE(&test_accum));
    master_test_suite().add(
        BOOST_PARAM_TEST_CASE(&test_gsl_rng, counts.begin(), counts.end()-1));
    master_test_suite().add(
        BOOST_PARAM_TEST_CASE(&test_rand48_gpu, counts.begin(), counts.end()));

    return 0;
}

static int _dummy = init_unit_test_suite();
