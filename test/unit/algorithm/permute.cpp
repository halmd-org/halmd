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

#include <halmd/config.hpp>

#define BOOST_TEST_MODULE permute
#include <boost/test/unit_test.hpp>

#include <halmd/algorithm/host/permute.hpp>
#include <halmd/numeric/accumulator.hpp>
#include <halmd/utility/scoped_timer.hpp>
#include <halmd/utility/timer.hpp>
#include <test/tools/ctest.hpp>
#include <test/tools/init.hpp>

#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include <algorithm>
#include <random>
#include <vector>

/**
 * Returns array of uniformly distributed unsigned integers.
 */
static std::vector<unsigned int> make_uniform_array(int count)
{
    std::vector<unsigned int> output(count);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<unsigned int> dist;
    std::generate(output.begin(), output.end(), [&]() {
        return dist(gen);
    });
    return std::move(output);
}

/**
 * Test generation of permutation using halmd::permute on host.
 */
static void test_permutation_host(int count, int repeat)
{
    std::vector<unsigned int> input_key(count);
    std::iota(input_key.begin(), input_key.end(), 0);
    {
        std::vector<unsigned int> const order = make_uniform_array(count);
        std::sort(
            input_key.begin()
          , input_key.end()
          , [&](unsigned int j, unsigned int k) {
                return order[j] < order[k];
          }
        );
    }
    std::vector<unsigned int> const& const_input_key = input_key;
    std::vector<unsigned int> const input_value = make_uniform_array(count);

    BOOST_TEST_MESSAGE( "  " << count << " elements" );
    BOOST_TEST_MESSAGE( "  " << repeat << " iterations" );

    halmd::accumulator<double> elapsed;
    for (int i = 0; i < repeat; ++i) {
        std::vector<unsigned int> output_value(
            input_value.begin()
          , input_value.end()
        );
        {
            halmd::scoped_timer<halmd::timer> t(elapsed);
            BOOST_CHECK( halmd::permute(
                output_value.begin()
              , output_value.end()
              , const_input_key.begin()) == const_input_key.end()
            );
        }
        auto key_to_value = [&](unsigned int key) {
            return input_value[key];
        };
        BOOST_CHECK_EQUAL_COLLECTIONS(
            boost::make_transform_iterator(input_key.begin(), key_to_value)
          , boost::make_transform_iterator(input_key.end(), key_to_value)
          , output_value.begin()
          , output_value.end()
        );
    }
    BOOST_TEST_MESSAGE( "  " << mean(elapsed) * 1e3 << " ± " << error_of_mean(elapsed) * 1e3 << " ms per iteration" );
}

HALMD_TEST_INIT( permute )
{
    using namespace boost::unit_test;

    for (int count : {1000000, 100000, 10000, 1024, 1000, 512, 100, 10, 2, 1}) {
        test_suite* ts = BOOST_TEST_SUITE( std::to_string(count) );
        framework::master_test_suite().add(ts);
        {
            int const repeat = std::max(100000 / count, 5);
            auto permutation_host = [=]() {
                test_permutation_host(count, repeat);
            };
            ts->add(BOOST_TEST_CASE( permutation_host ));
        }
    }
}
