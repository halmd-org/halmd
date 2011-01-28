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

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <halmd/utility/timer.hpp>

using namespace halmd;

//
// benchmark CUDA memory allocation
//
void performance(size_t size)
{
    const int count = 100;
    double alloc = 0;
    double copy = 0;

    try {
        cuda::vector<char> g_array1(size);

        halmd::timer timer;
        for (int i = 0; i < count; ++i) {
            // allocate memory
            timer.restart();
            cuda::vector<char> g_array2(size);
            alloc += timer.elapsed();

            // copy previous array
            timer.restart();
            g_array2 = g_array1;
            copy += timer.elapsed();
        }
        alloc /= count;
        copy /= count;
        BOOST_TEST_MESSAGE("  allocation of " << size << " bytes: " << alloc * 1e3 << " ms");
        BOOST_TEST_MESSAGE("  copying of " << size << " bytes: " << copy * 1e3 << " ms");
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
    for (size_t s = 1; s <= (1L << 32); s <<= 1) {
        sizes.push_back(s);
    }

    master_test_suite().add(BOOST_PARAM_TEST_CASE(&performance, sizes.begin(), sizes.end()));
}
