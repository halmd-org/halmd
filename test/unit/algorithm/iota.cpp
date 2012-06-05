/*
 * Copyright © 2012 Peter Colberg
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

#define BOOST_TEST_MODULE iota
#include <boost/test/unit_test.hpp>

#include <halmd/algorithm/gpu/iota.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>
#define HALMD_TEST_NO_LOGGING
#include <test/tools/ctest.hpp>
#include <test/tools/cuda.hpp>
#include <test/tools/init.hpp>

#include <boost/iterator/counting_iterator.hpp>

#include <algorithm>

/**
 * Test halmd::iota on GPU.
 */
template <typename size_type>
static void test_iota(
    size_type count
  , size_type value
  , size_type repeat
)
{
    BOOST_TEST_MESSAGE( "  " << count << " elements with first value " << value );
    BOOST_TEST_MESSAGE( "  " << repeat << " iterations" );

    halmd::accumulator<double> elapsed;
    for (size_type i = 0; i < repeat; ++i) {
        cuda::vector<size_type> g_output(count);
        cuda::memset(g_output.begin(), g_output.end(), 0);
        {
            halmd::scoped_timer<halmd::timer> t(elapsed);
            halmd::iota(g_output.begin(), g_output.end(), value);
        }
        cuda::host::vector<size_type> h_output(count);
        BOOST_CHECK( cuda::copy(
            g_output.begin()
          , g_output.end()
          , h_output.begin()) == h_output.end()
        );
        BOOST_CHECK_EQUAL_COLLECTIONS(
            h_output.begin()
          , h_output.end()
          , boost::make_counting_iterator(value)
          , boost::make_counting_iterator(value + count)
        );
    }
    BOOST_TEST_MESSAGE( "  " << mean(elapsed) * 1e3 << " ± " << error_of_mean(elapsed) * 1e3 << " ms per iteration" );
}

/**
 * Manually register parametrised test cases.
 */
HALMD_TEST_INIT( iota )
{
    using namespace boost::unit_test;

    for (unsigned int count : {1000000, 100000, 10000, 1024, 1000, 512, 100, 10, 2, 1}) {
        test_suite* ts = BOOST_TEST_SUITE( std::to_string(count) );
        framework::master_test_suite().add(ts);

        unsigned int const repeat = std::max(100 / count, 10u);

        for (unsigned int value : {0, 42}) {
            auto iota = [=]() {
                set_cuda_device device;
                test_iota(count, value, repeat);
            };
            ts->add(BOOST_TEST_CASE( iota ));
        }
    }
}
