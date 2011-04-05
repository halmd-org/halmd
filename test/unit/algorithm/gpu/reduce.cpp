/*
 * Copyright © 2010-2011  Peter Colberg and Felix Höfling
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

#define BOOST_TEST_MODULE reduce
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>

#include <boost/array.hpp>
#include <boost/iterator/counting_iterator.hpp>

#include <halmd/algorithm/gpu/reduce.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>
#include <halmd/utility/timer.hpp>
#include <test/tools/cuda.hpp>
#include <test/tools/init.hpp>

using namespace boost;
using namespace halmd;
using namespace halmd::algorithm::gpu;

BOOST_GLOBAL_FIXTURE( set_cuda_device );

/**
 * Test »32 bit integer arithmetic using double-single floating point (~44 bit).
 */
BOOST_AUTO_TEST_CASE( double_single_reduce )
{
    reduce<
        sum_            // reduce_transform
      , float           // input_type
      , float           // coalesced_input_type
      , dsfloat         // output_type
      , dsfloat         // coalesced_output_type
      , double          // host_output_type
    > sum;
    cuda::host::vector<float> h_v(make_counting_iterator(1), make_counting_iterator(12345679));
    cuda::vector<float> g_v(h_v.size());
    cuda::copy(h_v, g_v);
    BOOST_CHECK_EQUAL(int64_t(sum(g_v)), 12345678LL * 12345679 / 2);
}

/**
 * benchmark reduce kernel
 */
void performance(size_t size)
{
    const int count_local = 20;
    const int count_global = 100;
    const int blocks = 30;
    const int threads = 64 << 3;

    try {
        cuda::host::vector<float> h_v(make_counting_iterator(size_t(1)), make_counting_iterator(size));
        cuda::vector<float> g_v(h_v.size());
        cuda::copy(h_v, g_v);

        array<accumulator<double>, 2> elapsed;
        for (int i = 0; i < count_local; ++i) {
            halmd::timer timer;
            // initialise kernel and allocate internal memory
            reduce<
                sum_            // reduce_transform
              , float           // input_type
              , float           // coalesced_input_type
              , dsfloat         // output_type
              , dsfloat         // coalesced_output_type
              , double          // host_output_type
            > sum_local(blocks, threads);
            //
            double result = sum_local(g_v);
            elapsed[0](timer.elapsed());
            BOOST_CHECK_EQUAL(uint64_t(result), uint64_t(size - 1) * size / 2);
        }

        // pre-initialise kernel
        reduce<
            sum_              // reduce_transform
            , float           // input_type
            , float           // coalesced_input_type
            , dsfloat         // output_type
            , dsfloat         // coalesced_output_type
            , double          // host_output_type
        > sum_global(blocks, threads);
        for (int i = 0; i < count_global; ++i) {
            halmd::timer timer;
            double result = sum_global(g_v);
            elapsed[1](timer.elapsed());
            BOOST_CHECK_EQUAL(uint64_t(result), uint64_t(size - 1) * size / 2);
        }

        BOOST_TEST_MESSAGE("  summation of " << size << " floats: "
            << mean(elapsed[0]) * 1e3 << " ± " << error_of_mean(elapsed[0]) * 1e3 << " ms (local), "
            << mean(elapsed[1]) * 1e3 << " ± " << error_of_mean(elapsed[1]) * 1e3 << " ms (global)"
        );
    }
    catch (cuda::error const& e) {
        BOOST_CHECK(e.err == cudaErrorMemoryAllocation);
        BOOST_TEST_MESSAGE("Unsufficient device memory. Skip test for array of " << size << " bytes");
    }
    catch (std::bad_alloc const& e) {
        BOOST_TEST_MESSAGE("Unsufficient host memory. Skip test for array of " << size << " bytes");
    }
}

HALMD_TEST_INIT( init_unit_test_suite )
{
    using namespace boost::unit_test::framework;

    std::vector<size_t> sizes;
    for (size_t s = 1; s <= (1L << 29); s <<= 1) {
        sizes.push_back(s);
    }

    master_test_suite().add(BOOST_PARAM_TEST_CASE(&performance, sizes.begin(), sizes.end()));
}
