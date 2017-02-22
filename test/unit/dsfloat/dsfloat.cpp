/*
 * Copyright © 2016 Daniel Kirchner
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

#define BOOST_TEST_MODULE dsfloat
#include <boost/test/unit_test.hpp>

#include <halmd/numeric/accumulator.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>
#define HALMD_TEST_NO_LOGGING
#include <test/tools/ctest.hpp>
#include <test/tools/cuda.hpp>
#include <test/tools/init.hpp>
#include <test/unit/dsfloat/dsfloat.hpp>

/**
 * Test dsfloat
 */
static void test_dsfloat_performance() {
    BOOST_TEST_MESSAGE("dsfloat performance test");
    unsigned int memsize = 4096 * 1024;

    cuda::vector<float4> data(memsize);
    data.reserve(memsize * 2);
    cuda::memset(data.begin(), data.begin() + data.capacity(), 0);

    halmd::fixed_vector<halmd::dsfloat, 3> increment;
    increment[0] = increment[1] = increment[2] = 0.1;

    auto dim = cuda::config(memsize / 128, 128);

    int iterations = 256;

    {
        halmd::accumulator<double> elapsed;

        for (auto i = 0; i < iterations; i++) {
            cuda::configure(dim.grid, dim.block);
            {
                halmd::scoped_timer<halmd::timer> t(elapsed);
                dsfloat_kernel_wrapper::kernel.test1(data, increment);
                cuda::thread::synchronize();
            }
        }
        BOOST_TEST_MESSAGE("  " << mean(elapsed) * 1e3 << " ± " << error_of_mean(elapsed) * 1e3 << " ms per iteration");
    }
    {
        halmd::accumulator<double> elapsed;
        halmd::dsfloat_vector<float4> data(memsize);

        for (auto i = 0; i < iterations; i++) {
            cuda::configure(dim.grid, dim.block);
            {
                halmd::scoped_timer<halmd::timer> t(elapsed);
                dsfloat_kernel_wrapper::kernel.test2(data.data(), increment);
                cuda::thread::synchronize();
            }
        }
        BOOST_TEST_MESSAGE("  " << mean(elapsed) * 1e3 << " ± " << error_of_mean(elapsed) * 1e3 << " ms per iteration");
    }
}

static void test_dsfloat_overload()
{
    auto dim = cuda::config(256/128, 128);

    BOOST_TEST_MESSAGE( "dsfloat overload test" );
    {
        cuda::host::vector<halmd::fixed_vector<float,3>> result1 (256);
        cuda::host::vector<halmd::fixed_vector<float,3>> result2 (256);

        halmd::fixed_vector<float,3> float_increment (0.1f);
        halmd::fixed_vector<halmd::dsfloat,3> dsfloat_increment (0.1);

        cuda::vector<float4> data (256);
        data.reserve(256*2);

        cuda::memset(data.begin(), data.begin()+data.capacity(), 0);

        cuda::configure(dim.grid, dim.block);
        dsfloat_kernel_overloaded_wrapper<float>::kernel.test(data, float_increment);
        cuda::thread::synchronize();

        {
            cuda::host::vector<float4> tmp (256);
            cuda::copy(data.begin(), data.end(), tmp.begin());
            int ignored;
            for (auto i = 0; i < 256; i++) {
                halmd::tie(result1[i], ignored) <<= tmp[i];
            }
        }

        cuda::memset(data.begin(), data.begin()+data.capacity(), 0);

        halmd::dsfloat_vector<float4> dsdata (256);
        cuda::vector<float4> &dsdata_float4 = dsdata;

        cuda::memset(dsdata_float4.begin(), dsdata_float4.begin() + dsdata_float4.capacity(), 0);

        cuda::configure(dim.grid, dim.block);
        dsfloat_kernel_overloaded_wrapper<halmd::dsfloat>::kernel.test(dsdata, dsfloat_increment);
        cuda::thread::synchronize();

        {
            cuda::host::vector<float4> tmp (256);
            cuda::copy(dsdata_float4.begin(), dsdata_float4.end(), tmp.begin());
            int ignored;
            for (auto i = 0; i < 256; i++) {
                halmd::tie(result2[i], ignored) <<= tmp[i];
            }
        }

        using namespace halmd;
        BOOST_CHECK_EQUAL_COLLECTIONS(
             result1.begin(), result1.end()
           , result2.begin(), result2.end()
        );
    }

}

BOOST_AUTO_TEST_CASE( dsfloat )
{
    test_dsfloat_overload();
    test_dsfloat_performance();
}
