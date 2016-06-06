/*
 * Copyright © 2008  Peter Colberg
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

#define BOOST_TEST_MODULE rand48
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <vector>

#include <halmd/random/gpu/rand48.hpp>
#include <halmd/random/gpu/random_kernel.hpp>
#include <halmd/utility/timer.hpp>
#include <test/tools/ctest.hpp>
#include <test/tools/cuda.hpp>

//
// Parallel GPU rand48 random number generator test
//

BOOST_GLOBAL_FIXTURE( set_cuda_device );

BOOST_AUTO_TEST_CASE( compare_variates )
{
    const unsigned blocks = 64;
    const unsigned threads = 128;
    const unsigned seed = time(NULL);
    const unsigned count = 10000000u;

    BOOST_TEST_MESSAGE("number of integers: " << count);
    BOOST_TEST_MESSAGE("blocks: " << blocks);
    BOOST_TEST_MESSAGE("threads: " << threads);
    BOOST_TEST_MESSAGE("seed: " << seed);

    using halmd::random::gpu::rand48;

    // seed GPU random number generator
    rand48 rng(blocks, threads);
    halmd::timer timer;
    timer.restart();
    rng.seed(seed);
    double elapsed = timer.elapsed();

    BOOST_TEST_MESSAGE("seed GPU time: " << std::fixed << std::setprecision(3)
                       << elapsed * 1e3 << " ms");

    // parallel GPU rand48
    cuda::vector<uint> g_array(count);
    timer.restart();
    cuda::configure(rng.dim.grid, rng.dim.block);
    halmd::random::gpu::get_random_kernel<rand48::rng_type>().get(g_array, g_array.size(), rng.rng());
    cuda::thread::synchronize();
    elapsed = timer.elapsed();

    BOOST_TEST_MESSAGE("rand48 GPU time: " << std::fixed << std::setprecision(3)
                       << elapsed * 1e3 << " ms");

    cuda::host::vector<uint> h_array(count);
    cuda::copy(g_array, h_array);

    // serial GNU C library rand48
    std::vector<uint> h_array2(count);
    srand48(seed);
    timer.restart();
    std::generate(h_array2.begin(), h_array2.end(), mrand48);
    elapsed = timer.elapsed();

    BOOST_TEST_MESSAGE("rand48 CPU time: " << std::fixed << std::setprecision(3)
                       << elapsed * 1e3 << " ms");

    // compare GPU and CPU variates
    BOOST_CHECK_MESSAGE(std::equal(h_array.begin(), h_array.end(), h_array2.begin()),
                        "random variates mismatch");

    unsigned mismatches = 0;
    for (unsigned i=0; i < count; ++i) {
        if (h_array[i] != h_array2[i]) {
            ++mismatches;
            if  (mismatches < 10) {
                BOOST_TEST_MESSAGE(i << " " << h_array[i] << " " << h_array2[i]);
            }
        }
    }
    if (mismatches) {
        BOOST_TEST_MESSAGE(mismatches << " mismatches out of " << count);
    }

}
