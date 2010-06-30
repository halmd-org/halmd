/*
 * Copyright © 2010  Felix Höfling
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

#include <halmd/numeric/host/accumulator.hpp>
#include <halmd/random/host/gsl_rng.hpp>
#include <halmd/random/gpu/rand48.hpp>

//
// test random number generators in combination with accumulator class
//

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
        halmd::random::gpu::rand48 rng(cuda::config(BLOCKS, THREADS));
        rng.set(seed);
        cuda::thread::synchronize();

        // parallel GPU rand48
        cuda::vector<float> g_array(count);
        rng.uniform(g_array);
        cuda::thread::synchronize();

        cuda::host::vector<float> h_array(count);
        cuda::copy(g_array, h_array);

        halmd::accumulator<double> a;
        for (unsigned long i=0; i < count; i++) {
            a += h_array[i];
         }

        // check count, mean, and variance
        BOOST_CHECK_EQUAL(a.count(), count);

        // mean = 1/2, variance = 1/12, n-th moment is 1/(n+1)
        // use tolerance = 3 sigma, so the test passes with 99% probability
        double val = 0.5;
        double tol = 3 * a.std() / sqrt(count - 1);
        BOOST_CHECK_CLOSE_FRACTION(a.mean(), val, tol / val);

        // Var(ΣX^2/N) = E(X^4)/N
        val = 1./12;
        tol = 1 * sqrt(1. / (count - 1) * (1./5));
        BOOST_CHECK_CLOSE_FRACTION(a.var(), val, tol / val);

        // TODO: Kolmogorov-Smirnov test
        // see Knuth, vol. 2, ch. 3.3.B
        // lookup algorithm without sorting by Gonzalez et al.
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

    // seed host random number generator
    halmd::random::host::gfsr4 rng;
    rng.set(seed);

    // test uniform distribution
    halmd::accumulator<double> a;
    for (unsigned i=0; i < count; i++) {
        a += rng.uniform();
    }

    // check count, mean, and variance
    BOOST_CHECK_EQUAL(a.count(), count);

    // mean = 1/2, variance = 1/12, n-th moment is 1/(n+1)
    // use tolerance = 3 sigma, so the test passes with 99% probability
    double val = 0.5;
    double tol = 3 * a.std() / sqrt(count - 1);
    BOOST_CHECK_CLOSE_FRACTION(a.mean(), val, tol / val);

    // Var(ΣX^2/N) = E(X^4)/N
    val = 1./12;
    tol = 1 * sqrt(1. / (count - 1) * (1./5));
    BOOST_CHECK_CLOSE_FRACTION(a.var(), val, tol / val);

    // test Gaussian distribution
    BOOST_REQUIRE_MESSAGE(count % 2 == 0, "Box-Muller algorithm requires even number of samples");
    a.clear();
    halmd::accumulator<double> a3, a4;
    for (unsigned i=0; i < count; i+=2) {
        double x1, x2, tmp;
        rng.gaussian(x1, x2, 1.);
        a += x1;
        tmp = x1 * x1;
        a3 += x1 * tmp;
        a4 += tmp * tmp;

        a += x2;
        tmp = x2 * x2;
        a3 += x2 * tmp;
        a4 += tmp * tmp;
    }

    // mean = 0, std = 1
    BOOST_CHECK_EQUAL(a.count(), count);
    tol = 3 * a.std() / sqrt(count - 1);         // tolerance = 3 sigma (passes in 99% of all cases)
    BOOST_CHECK_SMALL(a.mean(), tol);
    val = 1;
    tol = 3 * sqrt( 1. / (count - 1) * 2);       // <X⁴> = 3 <X²>² = 3 ⇒ Var(X²) = 2
    BOOST_CHECK_CLOSE_FRACTION(a.var(), val, tol / val);

    // higher moments
    tol = 3 * a3.std() / sqrt(count - 1);        // <X³> = 0
    BOOST_CHECK_SMALL(a3.mean(), tol);
    val = 3;                                     // <X⁴> = 3
    tol = 3 * a4.std() / sqrt(count - 1);
    BOOST_CHECK_CLOSE_FRACTION(a4.mean(), val, tol / val);
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
        BOOST_PARAM_TEST_CASE(&test_gsl_rng, counts.begin(), counts.end()-2));
    master_test_suite().add(
        BOOST_PARAM_TEST_CASE(&test_rand48_gpu, counts.begin(), counts.end()));

    return 0;
}

static int _dummy = init_unit_test_suite();
