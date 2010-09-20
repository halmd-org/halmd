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

#define BOOST_TEST_MODULE test_random_distributions
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>

#include <boost/random/uniform_01.hpp>
#include <boost/random/normal_distribution.hpp>
#include <time.h>
#include <stdexcept>

#include <cuda_wrapper.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <halmd/random/gpu/random_kernel.hpp>
#include <halmd/random/host/gsl_rng.hpp>
#include <halmd/random/gpu/rand48.hpp>
#include <halmd/test/tools/cuda.hpp>

//
// test random number generators in combination with accumulator class
//

const unsigned BLOCKS = 64;
const unsigned THREADS = 128;

BOOST_GLOBAL_FIXTURE( set_cuda_device );

void test_rand48_gpu( unsigned long n )
{
    unsigned seed = time(NULL);
    using halmd::random::gpu::rand48;

    try {
        // seed GPU random number generator
        rand48 rng(BLOCKS, THREADS);
        rng.seed(seed);
        cuda::copy(rng.rng(), halmd::random::gpu::get_random_kernel<rand48::rng_type>().rng);

        // parallel GPU rand48
        cuda::vector<float> g_array(n);
        cuda::configure(rng.dim.grid, rng.dim.block);
        halmd::random::gpu::get_random_kernel<rand48::rng_type>().uniform(g_array, g_array.size());
        cuda::thread::synchronize();

        cuda::host::vector<float> h_array(n);
        cuda::copy(g_array, h_array);

        halmd::accumulator<double> a;
        for (unsigned long i=0; i < n; i++) {
            a(h_array[i]);
         }

        // check count, mean, and variance
        BOOST_CHECK_EQUAL(count(a), n);

        // mean = 1/2, variance = 1/12, n-th moment is 1/(n+1)
        // use tolerance = 3 sigma, so the test passes with 99% probability
        double val = 0.5;
        double tol = 3 * sigma(a) / sqrt(n - 1);
        BOOST_CHECK_CLOSE_FRACTION(mean(a), val, tol / val);

        // Var(ΣX^2/N) = E(X^4)/N
        val = 1./12;
        tol = 1 * sqrt(1. / (n - 1) * (1./5));
        BOOST_CHECK_CLOSE_FRACTION(variance(a), val, tol / val);

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

void test_gsl_rng( unsigned long n )
{
    typedef halmd::random::host::gfsr4 RandomNumberGenerator;

    unsigned seed = time(NULL);

    // seed host random number generator
    RandomNumberGenerator rng(seed);
    // uniform distribution in [0.0, 1.0)
    boost::uniform_01<RandomNumberGenerator&, double> uniform_01(rng);

    // test uniform distribution
    halmd::accumulator<double> a;
    for (unsigned i=0; i < n; i++) {
        a(uniform_01());
    }

    // check count, mean, and variance
    BOOST_CHECK_EQUAL(count(a), n);

    // mean = 1/2, variance = 1/12, n-th moment is 1/(n+1)
    // use tolerance = 3 sigma, so the test passes with 99% probability
    double val = 0.5;
    double tol = 3 * sigma(a) / sqrt(n - 1);
    BOOST_CHECK_CLOSE_FRACTION(mean(a), val, tol / val);

    // Var(ΣX^2/N) = E(X^4)/N
    val = 1./12;
    tol = 1 * sqrt(1. / (n - 1) * (1./5));
    BOOST_CHECK_CLOSE_FRACTION(variance(a), val, tol / val);

    // test Gaussian distribution
    boost::normal_distribution<double> normal;
    BOOST_CHECK_EQUAL(normal.mean(), 0.0);
    BOOST_CHECK_EQUAL(normal.sigma(), 1.0);
    a = halmd::accumulator<double>();
    halmd::accumulator<double> a3, a4;
    for (unsigned i=0; i < n; i++) {
        double x = normal(uniform_01);
        a(x);
        double x2 = x * x;
        a3(x * x2);
        a4(x2 * x2);
    }

    // mean = 0, std = 1
    BOOST_CHECK_EQUAL(count(a), n);
    tol = 3 * sigma(a) / sqrt(n - 1);         // tolerance = 3 sigma (passes in 99% of all cases)
    BOOST_CHECK_SMALL(mean(a), tol);
    val = 1;
    tol = 3 * sqrt( 1. / (n - 1) * 2);       // <X⁴> = 3 <X²>² = 3 ⇒ Var(X²) = 2
    BOOST_CHECK_CLOSE_FRACTION(variance(a), val, tol / val);

    // higher moments
    tol = 3 * sigma(a3) / sqrt(n - 1);        // <X³> = 0
    BOOST_CHECK_SMALL(mean(a3), tol);
    val = 3;                                     // <X⁴> = 3
    tol = 3 * sigma(a4) / sqrt(n - 1);
    BOOST_CHECK_CLOSE_FRACTION(mean(a4), val, tol / val);
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

    master_test_suite().add(
        BOOST_PARAM_TEST_CASE(&test_gsl_rng, counts.begin(), counts.end()-2));
    master_test_suite().add(
        BOOST_PARAM_TEST_CASE(&test_rand48_gpu, counts.begin(), counts.end()));

    return 0;
}

static int _dummy = init_unit_test_suite();
