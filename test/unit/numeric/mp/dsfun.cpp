/*
 * Copyright © 2009  Peter Colberg
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

#define BOOST_TEST_MODULE dsfun
#include <boost/test/unit_test.hpp>

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>
#include <test/tools/ctest.hpp>
#include <test/tools/cuda.hpp>

#define BLOCKS 4096
#define THREADS 128

using namespace halmd;

extern cuda::function<void (dsfloat const*, dsfloat const*, dsfloat*)> kernel_add;
extern cuda::function<void (dsfloat const*, dsfloat const*, dsfloat*)> kernel_sub;
extern cuda::function<void (dsfloat const*, dsfloat const*, dsfloat*)> kernel_mul;
extern cuda::function<void (float const*, float const*, dsfloat*)> kernel_mulss;
extern cuda::function<void (dsfloat const*, dsfloat const*, dsfloat*)> kernel_div;
extern cuda::function<void (dsfloat const*, dsfloat*)> kernel_sqrt;

BOOST_GLOBAL_FIXTURE( set_cuda_device )

/**
 * DSFUN addition test
 */
BOOST_AUTO_TEST_CASE(test_dsfun_add)
{
    cuda::config dim(BLOCKS, THREADS);
    cuda::host::vector<double> h_a(BLOCKS * THREADS);
    cuda::host::vector<double> h_b(h_a.size());
    cuda::host::vector<double> h_c(h_a.size());
    cuda::host::vector<dsfloat> h_dsp(h_a.size());
    cuda::vector<dsfloat> g_a(h_dsp.size());
    cuda::vector<dsfloat> g_b(h_dsp.size());
    cuda::vector<dsfloat> g_c(h_dsp.size());

    // generate double-precision random numbers on host
    srand48(42);
    std::generate(h_a.begin(), h_a.end(), drand48);
    // convert to double-single precision numbers
    std::copy(h_a.begin(), h_a.end(), h_dsp.begin());
    cuda::copy(h_dsp, g_a);

    std::generate(h_b.begin(), h_b.end(), drand48);
    std::copy(h_b.begin(), h_b.end(), h_dsp.begin());
    cuda::copy(h_dsp, g_b);

    cuda::memset(g_c, 0);
    cuda::configure(dim.grid, dim.block);
    kernel_add(g_a, g_b, g_c);
    cuda::thread::synchronize();

    cuda::copy(g_c, h_dsp);
    std::copy(h_dsp.begin(), h_dsp.end(), h_c.begin());

    for (size_t i = 0; i < h_c.size(); ++i) {
        BOOST_CHECK_CLOSE_FRACTION(h_c[i], h_a[i] + h_b[i], 1e-14);
    }
}

/**
 * DSFUN addition test
 */
BOOST_AUTO_TEST_CASE(test_dsfun_sub)
{
    cuda::config dim(BLOCKS, THREADS);
    cuda::host::vector<double> h_a(BLOCKS * THREADS);
    cuda::host::vector<double> h_b(h_a.size());
    cuda::host::vector<double> h_c(h_a.size());
    cuda::host::vector<dsfloat> h_dsp(h_a.size());
    cuda::vector<dsfloat> g_a(h_dsp.size());
    cuda::vector<dsfloat> g_b(h_dsp.size());
    cuda::vector<dsfloat> g_c(h_dsp.size());

    // generate double-precision random numbers on host
    srand48(42);
    std::generate(h_a.begin(), h_a.end(), drand48);
    // convert to double-single precision numbers
    std::copy(h_a.begin(), h_a.end(), h_dsp.begin());
    cuda::copy(h_dsp, g_a);

    std::generate(h_b.begin(), h_b.end(), drand48);
    std::copy(h_b.begin(), h_b.end(), h_dsp.begin());
    cuda::copy(h_dsp, g_b);

    cuda::memset(g_c, 0);
    cuda::configure(dim.grid, dim.block);
    kernel_sub(g_a, g_b, g_c);
    cuda::thread::synchronize();

    cuda::copy(g_c, h_dsp);
    std::copy(h_dsp.begin(), h_dsp.end(), h_c.begin());

    for (size_t i = 0; i < h_c.size(); ++i) {
        BOOST_CHECK_CLOSE_FRACTION(h_c[i], h_a[i] - h_b[i], 1e-15);
    }
}

/**
 * DSFUN multiplication test
 */
BOOST_AUTO_TEST_CASE(test_dsfun_mul)
{
    cuda::config dim(BLOCKS, THREADS);
    cuda::host::vector<double> h_a(BLOCKS * THREADS);
    cuda::host::vector<double> h_b(h_a.size());
    cuda::host::vector<double> h_c(h_a.size());
    cuda::host::vector<dsfloat> h_dsp(h_a.size());
    cuda::vector<dsfloat> g_a(h_dsp.size());
    cuda::vector<dsfloat> g_b(h_dsp.size());
    cuda::vector<dsfloat> g_c(h_dsp.size());

    // generate double-precision random numbers on host
    srand48(42);
    std::generate(h_a.begin(), h_a.end(), drand48);
    // convert to double-single precision numbers
    std::copy(h_a.begin(), h_a.end(), h_dsp.begin());
    cuda::copy(h_dsp, g_a);

    std::generate(h_b.begin(), h_b.end(), drand48);
    std::copy(h_b.begin(), h_b.end(), h_dsp.begin());
    cuda::copy(h_dsp, g_b);

    cuda::memset(g_c, 0);
    cuda::configure(dim.grid, dim.block);
    kernel_mul(g_a, g_b, g_c);
    cuda::thread::synchronize();

    cuda::copy(g_c, h_dsp);
    std::copy(h_dsp.begin(), h_dsp.end(), h_c.begin());

    for (size_t i = 0; i < h_c.size(); ++i) {
        BOOST_CHECK_CLOSE_FRACTION(h_c[i], h_a[i] * h_b[i], 2e-14);
    }
}

/**
 * DSFUN single-single multiplication test
 */
BOOST_AUTO_TEST_CASE(test_dsfun_mulss)
{
    cuda::config dim(BLOCKS, THREADS);
    cuda::host::vector<float> h_a(BLOCKS * THREADS);
    cuda::host::vector<float> h_b(h_a.size());
    cuda::host::vector<double> h_c(h_a.size());
    cuda::host::vector<dsfloat> h_dsp(h_a.size());
    cuda::vector<float> g_a(h_dsp.size());
    cuda::vector<float> g_b(h_dsp.size());
    cuda::vector<dsfloat> g_c(h_dsp.size());

    srand48(42);
    std::generate(h_a.begin(), h_a.end(), drand48);
    cuda::copy(h_a, g_a);
    std::generate(h_b.begin(), h_b.end(), drand48);
    cuda::copy(h_b, g_b);

    cuda::memset(g_c, 0);
    cuda::configure(dim.grid, dim.block);
    kernel_mulss(g_a, g_b, g_c);
    cuda::thread::synchronize();

    cuda::copy(g_c, h_dsp);
    std::copy(h_dsp.begin(), h_dsp.end(), h_c.begin());

    for (size_t i = 0; i < h_c.size(); ++i) {
        BOOST_CHECK_CLOSE_FRACTION(h_c[i], (double) h_a[i] * h_b[i], 1e-14);
    }
}

/**
 * DSFUN division test
 */
BOOST_AUTO_TEST_CASE(test_dsfun_div)
{
    cuda::config dim(BLOCKS, THREADS);
    cuda::host::vector<double> h_a(BLOCKS * THREADS);
    cuda::host::vector<double> h_b(h_a.size());
    cuda::host::vector<double> h_c(h_a.size());
    cuda::host::vector<dsfloat> h_dsp(h_a.size());
    cuda::vector<dsfloat> g_a(h_dsp.size());
    cuda::vector<dsfloat> g_b(h_dsp.size());
    cuda::vector<dsfloat> g_c(h_dsp.size());

    // generate double-precision random numbers on host
    srand48(42);
    std::generate(h_a.begin(), h_a.end(), drand48);
    // convert to double-single precision numbers
    std::copy(h_a.begin(), h_a.end(), h_dsp.begin());
    cuda::copy(h_dsp, g_a);

    std::generate(h_b.begin(), h_b.end(), drand48);
    std::copy(h_b.begin(), h_b.end(), h_dsp.begin());
    cuda::copy(h_dsp, g_b);

    cuda::memset(g_c, 0);
    cuda::configure(dim.grid, dim.block);
    kernel_div(g_a, g_b, g_c);
    cuda::thread::synchronize();

    cuda::copy(g_c, h_dsp);
    std::copy(h_dsp.begin(), h_dsp.end(), h_c.begin());

    for (size_t i = 0; i < h_c.size(); ++i) {
        BOOST_CHECK_CLOSE_FRACTION(h_c[i], h_a[i] / h_b[i], 5e-14);
    }
}

/**
 * DSFUN square-root test
 */
BOOST_AUTO_TEST_CASE(test_dsfun_sqrt)
{
    cuda::config dim(BLOCKS, THREADS);
    cuda::host::vector<double> h_a(BLOCKS * THREADS);
    cuda::host::vector<double> h_b(h_a.size());
    cuda::host::vector<dsfloat> h_dsp(h_a.size());
    cuda::vector<dsfloat> g_a(h_dsp.size());
    cuda::vector<dsfloat> g_b(h_dsp.size());

    // generate double-precision random numbers on host
    srand48(42);
    std::generate(h_a.begin(), h_a.end(), drand48);
    // convert to double-single precision numbers
    std::copy(h_a.begin(), h_a.end(), h_dsp.begin());
    cuda::copy(h_dsp, g_a);

    cuda::memset(g_b, 0);
    cuda::configure(dim.grid, dim.block);
    kernel_sqrt(g_a, g_b);
    cuda::thread::synchronize();

    cuda::copy(g_b, h_dsp);
    std::copy(h_dsp.begin(), h_dsp.end(), h_b.begin());

    for (size_t i = 0; i < h_b.size(); ++i) {
        BOOST_CHECK_CLOSE_FRACTION(h_b[i], std::sqrt(h_a[i]), 2e-13);
    }
}
