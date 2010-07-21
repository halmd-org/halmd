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
#define BOOST_TEST_MODULE test_radix_sort
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>

#include <algorithm>
#include <boost/array.hpp>
#include <boost/assign.hpp>
#include <cstdlib>
#include <deque>
#include <iomanip>
#include <vector>

#include <halmd/algorithm/gpu/radix_sort.hpp>
#include <halmd/random/gpu/rand48.hpp>
#include <halmd/random/gpu/random_kernel.hpp>
#include <halmd/test/tools/cuda.hpp>
#include <halmd/util/timer.hpp>

BOOST_GLOBAL_FIXTURE( set_cuda_device );

void compare_radix_sort( size_t count )
{
    const unsigned threads = 128;
    const unsigned seed = 42;

    halmd::high_resolution_timer start, stop;
    using halmd::random::gpu::rand48;

    // generate array of random integers in [0, 2^32-1] on GPU
    cuda::vector<uint> g_array(count);
    cuda::host::vector<uint> h_array(count);
    rand48 rng((count + threads - 1) / threads, threads);
    rng.seed(seed);
    cuda::copy(rng.rng(), halmd::random::gpu::get_random_kernel<rand48::rng_type>().rng);
    cuda::configure(rng.dim.grid, rng.dim.block);
    halmd::random::gpu::get_random_kernel<rand48::rng_type>().get(g_array, g_array.size());
    cuda::thread::synchronize();
    cuda::copy(g_array, h_array);

    // parallel radix sort
    halmd::algorithm::gpu::radix_sort<uint> sort(count, threads);
    cuda::vector<uint> g_dummy(count);
    cuda::host::vector<uint> h_array2(count);
    start.record();
    sort(g_array, g_dummy);
    cuda::thread::synchronize();
    stop.record();
    cuda::copy(g_array, h_array2);

    BOOST_TEST_MESSAGE("GPU time: " << std::fixed << std::setprecision(3)
                       << (stop - start) * 1e3 << " ms");

    using halmd::algorithm::gpu::BUCKET_SIZE;
    using halmd::algorithm::gpu::RADIX;

    // serial radix sort
    std::vector<uint> h_array3(count);
    boost::array<std::deque<uint>, BUCKET_SIZE> h_buckets;
    std::copy(h_array.begin(), h_array.end(), h_array3.begin());
    start.record();
    for (uint r = 0; r < 32; r += RADIX) {
        for (uint i = 0; i < count; ++i) {
            h_buckets[(h_array3[i] >> r) & 0xff].push_back(h_array3[i]);
        }
        h_array3.clear();
        for (uint b = 0; b < h_buckets.size(); ++b) {
            while (h_buckets[b].size()) {
                h_array3.push_back(h_buckets[b].front());
                h_buckets[b].pop_front();
            }
        }
    }
    stop.record();

    BOOST_TEST_MESSAGE("CPU time: " << std::fixed << std::setprecision(3)
                       << (stop - start) * 1e3 << " ms");

    // compare sorted result vectors
    BOOST_CHECK_MESSAGE(std::equal(h_array2.begin(), h_array2.end(), h_array3.begin()),
                        "GPU and CPU element mismatch");
}

int init_unit_test_suite()
{
    using namespace boost::unit_test::framework;
    using namespace boost::assign;

    std::vector<size_t> count;
    count += 1000000
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
        BOOST_PARAM_TEST_CASE(&compare_radix_sort, count.begin(), count.end()));

    return 0;
}

static int _dummy = init_unit_test_suite();
