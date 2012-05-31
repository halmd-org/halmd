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

#define BOOST_TEST_MODULE cache
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/reverse_iterator.hpp>
#include <vector>

#include <halmd/utility/cache.hpp>
#include <halmd/utility/raw_array.hpp>
#include <test/tools/ctest.hpp>

template <typename value_type, typename iterator_type>
static void test_cache(
    halmd::cache<value_type>& position
  , halmd::cache<>& position_cache
  , iterator_type const& first
  , iterator_type const& last
)
{
    BOOST_CHECK( position_cache == position );
    {
        halmd::cache_proxy<value_type const> position_proxy = position;
        BOOST_CHECK_EQUAL_COLLECTIONS(
            (*position_proxy).begin()
          , position_proxy->end()
          , first
          , last
        );
    }
    BOOST_CHECK( position == position_cache );
    {
        halmd::cache_proxy<value_type> position_proxy = position;
        BOOST_CHECK_EQUAL_COLLECTIONS(
            (*position_proxy).begin()
          , position_proxy->end()
          , first
          , last
        );
        std::reverse(position_proxy->begin(), position_proxy->end());
        BOOST_CHECK_EQUAL_COLLECTIONS(
            boost::make_reverse_iterator((*position_proxy).end())
          , boost::make_reverse_iterator(position_proxy->begin())
          , first
          , last
        );
    }
    BOOST_CHECK( position_cache != position );
    {
        halmd::cache_proxy<value_type const> position_proxy = position;
        BOOST_CHECK_EQUAL_COLLECTIONS(
            boost::make_reverse_iterator((*position_proxy).end())
          , boost::make_reverse_iterator(position_proxy->begin())
          , first
          , last
        );
    }
    BOOST_CHECK( position != position_cache );
}

/**
 * Test read and write access to cached value.
 */
BOOST_AUTO_TEST_CASE( read_write )
{
    // fill raw_array with ascending sequence of numbers
    std::vector<double> const input(
        boost::make_counting_iterator(0)
      , boost::make_counting_iterator(42)
    );
    // validates forwarding of parameters by the cache constructor
    halmd::cache<halmd::raw_array<double>> position(input.size());
    std::copy(
        input.begin()
      , input.end()
      , halmd::cache_proxy<halmd::raw_array<double>>(position)->begin()
    );

    // test read and write access, write reverses raw_array
    halmd::cache<> position_cache = position;
    test_cache(
        position
      , position_cache
      , input.begin()
      , input.end()
    );

    // test read and write access with reversed raw_array
    position_cache = position;
    test_cache(
        position
      , position_cache
      , input.rbegin()
      , input.rend()
    );
}

/**
 * Test cache observer with regard to cache lifetime.
 */
BOOST_AUTO_TEST_CASE( observer )
{
    halmd::cache<> position_cache;
    halmd::cache<double> position(M_PI);
    BOOST_CHECK( position_cache != position );

    position_cache = position;
    BOOST_CHECK( position_cache == position );

    BOOST_CHECK_EQUAL( *halmd::cache_proxy<double const>(position), M_PI );
    BOOST_CHECK( position == position_cache );

    position_cache = halmd::cache<>();
    BOOST_CHECK( position != position_cache );
}
