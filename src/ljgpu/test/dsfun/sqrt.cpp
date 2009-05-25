/* DSFUN square-root test
 *
 * Copyright (C) 2009  Peter Colberg
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

#include "sqrt_kernel.hpp"

#define BOOST_TEST_MODULE test_dsfun_sqrt
#include <boost/test/included/unit_test.hpp>

#define BLOCKS 4096
#define THREADS 128

BOOST_AUTO_TEST_CASE(test_dsfun_sqrt)
{
    cuda::config dim(BLOCKS, THREADS);
    cuda::host::vector<double> h_a(BLOCKS * THREADS);
    cuda::host::vector<double> h_b(h_a.size());
    cuda::host::vector<dfloat> h_dsp(h_a.size());
    cuda::vector<dfloat> g_a(h_dsp.size());
    cuda::vector<dfloat> g_b(h_dsp.size());

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
	BOOST_CHECK_CLOSE_FRACTION(h_b[i], std::sqrt(h_a[i]), 3e-14);
    }
}
