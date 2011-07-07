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

#define BOOST_TEST_MODULE fill_vs_memset
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>

#include <algorithm>
#include <string> // std::memset
#include <vector>

#include <halmd/utility/timer.hpp>
#include <test/tools/init.hpp>
#ifdef WITH_CUDA
# include <test/tools/cuda.hpp>
# include <test/performance/fill_kernel.hpp>
#endif

/**
 * This test compares zeroing the entries of std::vector< float >
 * via std::fill or via memset. The latter method requires IEEE 754 conformity.
 */
void host(size_t size)
{
    typedef float float_type;
    std::vector<float_type> a(size);
    size_t iterations = std::min(16384UL, (1UL << 28) / size);
    BOOST_TEST_MESSAGE("Measuring execution times over " << iterations << " iterations"
        << " for zeroing " << size << " floats"
    );

    double fill = 0;
    double memset = 0;
    for (size_t i = 0; i < iterations; ++i) {
        halmd::timer timer;
        std::fill(a.begin(), a.end(), 0);
        fill += timer.elapsed();

        timer.restart();
        std::memset(&a.front(), 0, sizeof(float_type) * a.size());
        memset += timer.elapsed();
    }
    BOOST_CHECK(a[0] == 0.f); // ensure IEEE 754 conformity

    fill /= iterations;
    memset /= iterations;
    BOOST_TEST_MESSAGE("std::fill: " << fill * 1e3 << " ms");
    BOOST_TEST_MESSAGE("std::memset: " << memset * 1e3 << " ms");
    double gain = fill / memset;
    BOOST_CHECK(gain > 1); // we expect memset to be faster
    BOOST_TEST_MESSAGE("gain factor of using memset: " << gain);
}

#ifdef WITH_CUDA

BOOST_GLOBAL_FIXTURE( set_cuda_device );

/**
 * This test compares zeroing the entries of cuda::vector< float >
 * via cuda::memset or via a CUDA kernel call.
 */
void gpu(size_t size)
{
    typedef float float_type;
    cuda::vector<float_type> a(size);
    size_t iterations = std::min(16384UL, (1UL << 28) / size);
    BOOST_TEST_MESSAGE("Measuring execution times over " << iterations << " iterations"
        << " for zeroing " << size << " floats"
    );

    unsigned int threads = 512;
    cuda::config dim((size + threads - 1) / threads, threads);
    if (dim.blocks_per_grid() > 65535) {
        return;
    }
    BOOST_TEST_MESSAGE("Using " << dim.threads_per_block() << " threads in each of " << dim.blocks_per_grid() << " blocks");

    double fill_loop = 0;
    double fill_if = 0;
    double memset = 0;

    for (size_t i = 0; i < iterations; ++i) {
        halmd::timer timer;
        cuda::configure(dim.grid, dim.block);
        fill_loop_kernel(a, a.size(), 0);
        cuda::thread::synchronize();
        fill_loop += timer.elapsed();

        timer.restart();
        cuda::configure(dim.grid, dim.block);
        fill_if_kernel(a, a.size(), 0);
        cuda::thread::synchronize();
        fill_if += timer.elapsed();

        timer.restart();
        cuda::memset(a, 0);
        memset += timer.elapsed();
    }

    fill_loop /= iterations;
    fill_if /= iterations;
    memset /= iterations;
    BOOST_TEST_MESSAGE("fill kernel (with loop): " << fill_loop * 1e3 << " ms");
    BOOST_TEST_MESSAGE("fill kernel (with if): " << fill_if * 1e3 << " ms");
    BOOST_TEST_MESSAGE("cuda::memset: " << memset * 1e3 << " ms");
    double gain = fill_if / memset;
    BOOST_CHECK(gain > 1); // we expect memset to be faster than fill_if
    BOOST_TEST_MESSAGE("gain factor of using memset over fill kernel (with if): " << gain);
}

#endif // WITH_CUDA

HALMD_TEST_INIT( init_unit_test_suite )
{
    using namespace boost::unit_test::framework;

    std::vector<size_t> sizes;
    for (size_t s = 512; s <= (1UL << 24); s <<= 1) {
        sizes.push_back(s);
    }

    master_test_suite().add(BOOST_PARAM_TEST_CASE(&host, sizes.begin(), sizes.end()));
#ifdef WITH_CUDA
    master_test_suite().add(BOOST_PARAM_TEST_CASE(&gpu, sizes.begin(), sizes.end()));
#endif
}
