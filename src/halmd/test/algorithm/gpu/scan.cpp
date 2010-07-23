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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_scan
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>

#include <algorithm>
#include <boost/assign.hpp>
#include <iomanip>
#include <vector>

#include <halmd/algorithm/gpu/scan.hpp>
#include <halmd/test/tools/cuda.hpp>
#include <halmd/utility/timer.hpp>

BOOST_GLOBAL_FIXTURE( set_cuda_device );

void compare_scan( size_t count )
{
    const unsigned threads = 256;

    BOOST_TEST_MESSAGE( "Scanning " << count << " elements");

    // generate array of ascending integers
    cuda::host::vector<uint> h_array(count);
    cuda::vector<uint> g_array(count);
    for (uint i = 0; i < count; ++i) {
        h_array[i] = i + 1;
    }
    cuda::copy(h_array, g_array);

    // parallel exclusive prefix sum
    halmd::algorithm::gpu::scan<uint> scan(count, threads);
    cuda::host::vector<uint> h_array2(count);
    halmd::utility::timer timer;
    scan(g_array);
    cuda::thread::synchronize();
    double elapsed = timer.elapsed();
    cuda::copy(g_array, h_array2);

    BOOST_TEST_MESSAGE("GPU time: " << std::fixed << std::setprecision(3)
                       << elapsed * 1e3 << " ms");

    // serial prefix sum
    std::vector<uint> h_array3(count);
    timer.restart();
    h_array3[0] = 0;
    for (uint i = 1; i < count; ++i) {
        h_array3[i] = h_array[i - 1] + h_array3[i - 1];
    }
    elapsed = timer.elapsed();

    BOOST_TEST_MESSAGE("CPU time: " << std::fixed << std::setprecision(3)
                       << elapsed * 1e3 << " ms");

    // compare results
    BOOST_CHECK_MESSAGE(std::equal(h_array2.begin(), h_array2.end(), h_array3.begin()),
                        "GPU and CPU prefix sum mismatch");
}

int init_unit_test_suite()
{
    using namespace boost::unit_test::framework;
    using namespace boost::assign;

    std::vector<size_t> count;
    count += 10000000
           , 1000000
           , 100000
           , 10000
           , 1024
           , 1000
           , 512
           , 100
           , 10
           , 2
           ;

    master_test_suite().add(
        BOOST_PARAM_TEST_CASE(&compare_scan, count.begin(), count.end()));

    return 0;
}

static int _dummy = init_unit_test_suite();
