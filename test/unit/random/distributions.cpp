/*
 * Copyright © 2010  Felix Höfling
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <halmd/config.hpp>

#define BOOST_TEST_MODULE distributions
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>

#include <cmath>
#include <ctime>
#include <stdexcept>

#include <halmd/numeric/accumulator.hpp>
#include <halmd/random/host/random.hpp>
#ifdef HALMD_WITH_GPU
# include <cuda_wrapper/cuda_wrapper.hpp>
# include <halmd/random/gpu/random_kernel.hpp>
# include <halmd/random/gpu/rand48.hpp>
# include <test/tools/cuda.hpp>
#endif
#include <test/tools/ctest.hpp>
#include <test/tools/init.hpp>

//
// test random number generators in combination with accumulator class
//

#ifdef HALMD_WITH_GPU

const unsigned BLOCKS = 64;
const unsigned THREADS = 128;

struct stochastic_test
{
    stochastic_test()
    {
        BOOST_TEST_MESSAGE("This test has a stochastic outcome. It shall fail less than once in 10000 passes.");
    }
};

BOOST_GLOBAL_FIXTURE( set_cuda_device );
BOOST_GLOBAL_FIXTURE( stochastic_test );

void test_rand48_gpu( unsigned long n )
{
    unsigned seed = time(NULL);
    using halmd::random::gpu::rand48;

    BOOST_TEST_MESSAGE("generate " << n << " uniformly distributed random numbers on the GPU");

    try {
        // seed GPU random number generator
        rand48 rng(BLOCKS, THREADS);
        rng.seed(seed);

        // parallel GPU rand48
        cuda::memory::device::vector<float> g_array(n);
        halmd::random::gpu::get_random_kernel<rand48::rng_type>().uniform.configure(
            rng.dim.grid, rng.dim.block);
        halmd::random::gpu::get_random_kernel<rand48::rng_type>().uniform(
            g_array, g_array.size(), rng.rng());
        cuda::thread::synchronize();

        cuda::memory::host::vector<float> h_array(n);
        cuda::copy(g_array.begin(), g_array.end(), h_array.begin());

        halmd::accumulator<double> a;
        for (unsigned long i=0; i < n; i++) {
            a(h_array[i]);
         }

        // check count, mean, and variance
        BOOST_CHECK_EQUAL(count(a), n);

        // mean = 1/2, variance = 1/12, N-th moment is 1/(N+1)
        // use tolerance = 4.5 sigma, so the test passes with 99.999% probability
        // (assuming a normal distribution of the measured value, which is true for large n)
        double val = 0.5;
        double tol = 4.5 * sigma(a) / std::sqrt(n - 1.);
        BOOST_CHECK_CLOSE_FRACTION(mean(a), val, tol / val);

        // Var(ΣX^2/N) = E(X^4)/N
        val = 1./12;
        tol = 1 * std::sqrt(1. / (n - 1) * (1./5));
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

#endif /* HALMD_WITH_GPU */

void test_host_random( unsigned long n )
{
    typedef halmd::random::host::random RandomNumberGenerator;

    unsigned seed = time(NULL);

    // seed host random number generator
    RandomNumberGenerator rng(seed);

    // test uniform distribution
    BOOST_TEST_MESSAGE("generate " << n << " uniformly distributed random numbers on the host");

    halmd::accumulator<double> a;
    for (unsigned i=0; i < n; i++) {
        a(rng.uniform<double>());
    }

    // check count, mean, and variance
    BOOST_CHECK_EQUAL(count(a), n);

    // mean = 1/2, variance = 1/12, N-th moment is 1/(N+1)
    // use tolerance = 4.5 sigma, so the test passes with 99.999% probability
    // (assuming a normal distribution of the measured value, which is true for large n)
    double val = 0.5;
    double tol = 4.5 * sigma(a) / std::sqrt(n - 1.);
    BOOST_CHECK_CLOSE_FRACTION(mean(a), val, tol / val);

    // Var(ΣX^2/N) = E(X^4)/N
    val = 1./12;
    tol = 1 * std::sqrt(1. / (n - 1) * (1./5));
    BOOST_CHECK_CLOSE_FRACTION(variance(a), val, tol / val);

    // test Gaussian distribution
    BOOST_TEST_MESSAGE("generate " << n << " normally distributed random numbers on the host");

    a = halmd::accumulator<double>();
    halmd::accumulator<double> a3, a4;
    for (unsigned i=0; i < n; i++) {
        // FIXME y is unused due to impractical host::random::normal design
        double x, y;
        std::tie(x, y) = rng.normal(1.0);
        a(x);
        double x2 = x * x;
        a3(x * x2);
        a4(x2 * x2);
    }

    // mean = 0, std = 1
    BOOST_CHECK_EQUAL(count(a), n);
    tol = 4.5 * sigma(a) / std::sqrt(n - 1.);     // tolerance = 4.5 sigma (passes in 99.999% of all cases)
    BOOST_CHECK_SMALL(mean(a), tol);
    val = 1;
    tol = 4.5 * std::sqrt( 1. / (n - 1) * 2);     // <X⁴> = 3 <X²>² = 3 ⇒ Var(X²) = 2
    BOOST_CHECK_CLOSE_FRACTION(variance(a), val, tol / val);

    // higher moments
    tol = 4.5 * sigma(a3) / std::sqrt(n - 1.);    // <X³> = 0
    BOOST_CHECK_SMALL(mean(a3), tol);
    val = 3;                                     // <X⁴> = 3
    tol = 4.5 * sigma(a4) / std::sqrt(n - 1.);
    BOOST_CHECK_CLOSE_FRACTION(mean(a4), val, tol / val);
}

HALMD_TEST_INIT( init_unit_test_suite )
{
    using namespace boost::unit_test::framework;

    std::vector<unsigned long> counts;
    counts.push_back(10000);
    counts.push_back(100000);
    counts.push_back(1000000);
    counts.push_back(10000000);
    counts.push_back(100000000);

    master_test_suite().add(
        BOOST_PARAM_TEST_CASE(&test_host_random, counts.begin(), counts.end()-2));
#ifdef HALMD_WITH_GPU
    master_test_suite().add(
        BOOST_PARAM_TEST_CASE(&test_rand48_gpu, counts.begin(), counts.end()));
#endif
}
