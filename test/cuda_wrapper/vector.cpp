/*
 * Copyright © 2011  Felix Höfling
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

#define BOOST_TEST_MODULE test_cuda_vector
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include <cmath>

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <halmd/utility/timer.hpp>

using namespace halmd;

/**
 * benchmark CUDA memory allocation
 *
 * Memory allocation on the device and of page-locked host memory is
 * tested and benchmarked. It turns out that calls to cudaMallocHost
 * have an overhead of several milliseconds.
 *
 * Caution: cudaMallocHost appears to prefer freezing your machine
 * instead of throwing an exception. The test allocates at most
 * 256 MB of memory.
 */
void performance(size_t size)
{
    const int count = std::max(8, static_cast<int>(256 / log2(size)));
    double alloc_device = 0;
    double alloc_host = 0;
    double copy = 0;

    try {
        halmd::timer timer;
        for (int i = 0; i < count; ++i) {
            // allocate device memory
            timer.restart();
            cuda::vector<char> g_array(size);
            alloc_device += timer.elapsed();

            // allocate page-locked host memory
            timer.restart();
            cuda::host::vector<char> h_array(size);
            alloc_host += timer.elapsed();

            // copy from device to host
            timer.restart();
            cuda::copy(g_array, h_array);
            copy += timer.elapsed();
        }
        alloc_device /= count;
        alloc_host /= count;
        copy /= count;
        BOOST_TEST_MESSAGE("  allocation of " << size << " bytes on the device: " << alloc_device * 1e3 << " ms");
        BOOST_TEST_MESSAGE("  allocation of " << size << " bytes on the host: " << alloc_host * 1e3 << " ms");
        BOOST_TEST_MESSAGE("  copying of " << size << " bytes (device → host): " << copy * 1e3 << " ms");
    }
    catch (cuda::error const& e) {
        BOOST_CHECK(e.err == cudaErrorMemoryAllocation);
        BOOST_TEST_MESSAGE("  skip test for chunk of " << size << " bytes");
    }
}

static void __attribute__((constructor)) init_unit_test_suite()
{
    using namespace boost::unit_test::framework;

    std::vector<size_t> sizes;
    for (size_t s = 4; s <= (1L << 28); s <<= 1) {
        sizes.push_back(s);
    }

    master_test_suite().add(BOOST_PARAM_TEST_CASE(&performance, sizes.begin(), sizes.end()));
}
