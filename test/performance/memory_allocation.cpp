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

#define BOOST_TEST_MODULE memory_allocation
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include <cmath>

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <halmd/utility/raw_array.hpp>
#include <halmd/utility/timer.hpp>
#include <test/tools/ctest.hpp>
#include <test/tools/cuda.hpp>
#include <test/tools/init.hpp>

BOOST_GLOBAL_FIXTURE( set_cuda_device );

using namespace halmd;

/**
 * benchmark CUDA and host memory allocation via 'vector' containers
 *
 * Memory allocation on the device and of page-locked host memory is
 * tested and benchmarked. It turns out that calls to cudaMallocHost
 * have an overhead of several milliseconds.
 *
 * Caution: cudaMallocHost appears to prefer freezing your machine
 * instead of throwing an exception. The test allocates at most
 * 256 MB of memory.
 */
void vector(size_t size)
{
    const int count = std::max(8, static_cast<int>(256 / log2(size)));
    double alloc_device = 0;
    double alloc_host = 0;
    double alloc_host_stl = 0;
    double alloc_host_raw = 0;
    double copy = 0;
    double copy_stl = 0;

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

            // allocate conventional host memory
            timer.restart();
            std::vector<char> stl_array(size);
            alloc_host_stl += timer.elapsed();

            // allocate conventional host memory without initialisation of elements
            timer.restart();
            halmd::raw_array<char> raw_array(size);
            alloc_host_raw += timer.elapsed();

            // copy from device to host
            timer.restart();
            cuda::copy(g_array, h_array);
            copy += timer.elapsed();

            // copy from host to host
            timer.restart();
            cuda::copy(h_array, stl_array);
            copy_stl += timer.elapsed();
        }
        alloc_device /= count;
        alloc_host /= count;
        alloc_host_stl /= count;
        alloc_host_raw /= count;
        copy /= count;
        copy_stl /= count;
        BOOST_TEST_MESSAGE("  allocation of " << size << " bytes on the device: " << alloc_device * 1e3 << " ms");
        BOOST_TEST_MESSAGE("  allocation of " << size << " bytes on the host (page-locked): " << alloc_host * 1e3 << " ms");
        BOOST_TEST_MESSAGE("  allocation of " << size << " bytes on the host (STL): " << alloc_host_stl * 1e3 << " ms");
        BOOST_TEST_MESSAGE("  allocation of " << size << " bytes on the host (raw): " << alloc_host_raw * 1e3 << " ms");
        BOOST_TEST_MESSAGE("  copying of " << size << " bytes (device → host): " << copy * 1e3 << " ms");
        BOOST_TEST_MESSAGE("  copying of " << size << " bytes (host → host): " << copy_stl * 1e3 << " ms");
    }
    catch (cuda::error const& e) {
        BOOST_CHECK(e.err == cudaErrorMemoryAllocation);
        BOOST_TEST_MESSAGE("  skip test for chunk of " << size << " bytes");
    }
}

HALMD_TEST_INIT( init_unit_test_suite )
{
    using namespace boost::unit_test::framework;

    std::vector<size_t> sizes;
    for (size_t s = 4; s <= (1L << 28); s <<= 1) {
        sizes.push_back(s);
    }

    master_test_suite().add(BOOST_PARAM_TEST_CASE(&vector, sizes.begin(), sizes.end()));
}
