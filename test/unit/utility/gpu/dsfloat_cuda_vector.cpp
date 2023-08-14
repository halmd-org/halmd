/*
 * Copyright © 2016 Daniel Kirchner
 * Copyright © 2017 Felix Höfling
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

#define BOOST_TEST_MODULE dsfloat_cuda_vector
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <limits>
#include <vector>

#include <halmd/numeric/accumulator.hpp>
#include <halmd/utility/gpu/dsfloat_cuda_vector.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>

#define HALMD_TEST_NO_LOGGING
#include <test/tools/ctest.hpp>
#include <test/tools/cuda.hpp>
#include <test/tools/init.hpp>
#include <test/unit/utility/gpu/dsfloat_cuda_vector_kernel.hpp>

using namespace halmd;

/**
 * Test performance of implementations: dsfloat_cuda_vector vs. float4 pointers
 */
BOOST_FIXTURE_TEST_CASE( performance, set_cuda_device )
{
    unsigned int memsize = 1024 * 1024;

    cuda::memory::device::vector<float4> data(memsize);
    data.reserve(memsize * 2);
    cuda::memset(data.begin(), data.begin() + data.capacity(), 0);

    fixed_vector<dsfloat, 3> increment;
    increment[0] = increment[1] = increment[2] = 0.1;

    auto dim = cuda::config(memsize / 128, 128);

    unsigned int iterations = 100;

    double mean_runtime_float4_ptr;
    double mean_runtime_dsfloat_ptr;
    {
        accumulator<double> elapsed;

        for (unsigned int i = 0; i < iterations; i++) {
            dsfloat_kernel_wrapper::kernel.test_float4_ptr.configure(dim.grid,
                dim.block);
            {
                scoped_timer<timer> t(elapsed);
                dsfloat_kernel_wrapper::kernel.test_float4_ptr(data, increment);
                cuda::thread::synchronize();
            }
        }
        mean_runtime_float4_ptr = mean(elapsed);
        BOOST_TEST_MESSAGE("  " << mean_runtime_float4_ptr * 1e3 << " ± " << error_of_mean(elapsed) * 1e3 << " ms per iteration");
    }
    {
        accumulator<double> elapsed;
        dsfloat_cuda_vector<float4> data(memsize);

        for (unsigned int i = 0; i < iterations; i++) {
            dsfloat_kernel_wrapper::kernel.test_dsfloat_ptr.configure(dim.grid,
                dim.block);
            {
                scoped_timer<timer> t(elapsed);
                dsfloat_kernel_wrapper::kernel.test_dsfloat_ptr(data.data(), increment);
                cuda::thread::synchronize();
            }
        }
        mean_runtime_dsfloat_ptr = mean(elapsed);
        BOOST_TEST_MESSAGE("  " << mean_runtime_dsfloat_ptr * 1e3 << " ± " << error_of_mean(elapsed) * 1e3 << " ms per iteration");
    }
    BOOST_CHECK_LE(mean_runtime_dsfloat_ptr, mean_runtime_float4_ptr * 1.05);
}

BOOST_FIXTURE_TEST_CASE( overloaded_kernel, set_cuda_device )
{
    // test parameters
    unsigned int const memsize = 256;
    auto dim = cuda::config(memsize / 128, 128);

    fixed_vector<float, 3> increment(std::numeric_limits<float>::epsilon());

    // result arrays
    std::vector<fixed_vector<float, 3>> result_float(memsize);
    std::vector<fixed_vector<float, 3>> result_dsfloat(memsize);

    // allocate GPU memory for 2×memsize float4, initialise to 0
    dsfloat_cuda_vector<float4> g_data(memsize);
    cuda::memory::device::vector<float4>& g_data_float4 = g_data;
    cuda::memset(g_data_float4.begin(), g_data_float4.begin() + g_data_float4.capacity(), 0);

    // assign 1's to first float4 array
    cuda::memory::host::vector<float4> h_data(memsize);
    std::fill(h_data.begin(), h_data.end(), float4{ 1, 1, 1, 1 });
    cuda::copy(h_data.begin(), h_data.end(), g_data_float4.begin());

    // use only float4 part of high significance
    dsfloat_kernel_overloaded_wrapper<float>::kernel.overloaded_test.configure(
        dim.grid, dim.block);
    dsfloat_kernel_overloaded_wrapper<float>::kernel.overloaded_test(
        g_data_float4, increment);
    cuda::thread::synchronize();

    // convert float4 to fixed_vector, subtract 1 to be sensitive to the last digits
    cuda::copy(g_data_float4.begin(), g_data_float4.end(), h_data.begin());
    for (size_t i = 0; i < h_data.size(); i++) {
        int ignored;
        tie(result_float[i], ignored) <<= h_data[i];
        result_float[i] -= fixed_vector<float, 3>(1);
    }

    // reset input array
    cuda::memset(g_data_float4.begin(), g_data_float4.begin() + g_data_float4.capacity(), 0);
    std::fill(h_data.begin(), h_data.end(), float4{ 1, 1, 1, 1 });
    cuda::copy(h_data.begin(), h_data.end(), g_data_float4.begin());

    // use full dsfloat data (high and low significant float4 arrays),
    // add a fraction of the previous increment many times so that the result should be the same,
    // but the high significant float4's are important
    unsigned int N = 1 << 8;
    for (unsigned int i = 0; i < N; ++i) {
        dsfloat_kernel_overloaded_wrapper<dsfloat>::kernel.overloaded_test.configure(
            dim.grid, dim.block);
        dsfloat_kernel_overloaded_wrapper<dsfloat>::kernel.overloaded_test(
            g_data, increment / N);
        cuda::thread::synchronize();
    }

    // convert float4 to fixed_vector, subtract 1 to be sensitive to the last digits
    cuda::copy(g_data_float4.begin(), g_data_float4.end(), h_data.begin());
    for (size_t i = 0; i < h_data.size(); i++) {
        int ignored;
        tie(result_dsfloat[i], ignored) <<= h_data[i];
        result_dsfloat[i] -= fixed_vector<float, 3>(1);
    }

    // compare results from the two routes
    BOOST_CHECK_EQUAL_COLLECTIONS(
         result_float.begin(), result_float.end()
       , result_dsfloat.begin(), result_dsfloat.end()
    );
}
