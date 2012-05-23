/*
 * Copyright © 2011-2012 Peter Colberg
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

#define BOOST_TEST_MODULE clock
#include <boost/test/unit_test.hpp>

#include <limits>

#include <halmd/mdsim/clock.hpp>
#include <test/tools/ctest.hpp>

/**
 * Construct instance of simulation clock.
 */
struct clock_fixture
{
    typedef halmd::mdsim::clock clock_type;
    typedef clock_type::step_type step_type;
    typedef clock_type::time_type time_type;

    /**
     * Returns relative tolerance for comparison of time values.
     */
    static time_type epsilon()
    {
       return std::numeric_limits<time_type>::epsilon();
    }

    clock_type clock;
    clock_type const& clock_const;

    clock_fixture() : clock_const(clock) {}
};

BOOST_FIXTURE_TEST_CASE( set_timestep, clock_fixture )
{
    BOOST_CHECK_EQUAL( clock_const.step(), 0u );
    BOOST_CHECK_CLOSE_FRACTION( clock_const.time(), 0., epsilon() );

    /** counter for number of calls to set_timestep slot */
    int calls = 0;
    /** value of timestep passed to set_timestep slot */
    time_type timestep = 0;

    clock.on_set_timestep([&](time_type value) {
        ++calls;
        timestep = value;
    });

    BOOST_CHECK_EQUAL( calls, 0 );
    BOOST_CHECK_CLOSE_FRACTION( timestep, 0., epsilon() );
    BOOST_CHECK_THROW( clock_const.timestep(), std::logic_error );

    BOOST_CHECK_NO_THROW( clock.set_timestep(0.001) );
    BOOST_CHECK_EQUAL( calls, 1 );
    BOOST_CHECK_CLOSE_FRACTION( timestep, 0.001, epsilon() );
    BOOST_CHECK_CLOSE_FRACTION( clock_const.timestep(), 0.001, epsilon() );

    BOOST_CHECK_NO_THROW( clock.set_timestep(0.002) );
    BOOST_CHECK_EQUAL( calls, 2 );
    BOOST_CHECK_CLOSE_FRACTION( timestep, 0.002, epsilon() );
    BOOST_CHECK_CLOSE_FRACTION( clock_const.timestep(), 0.002, epsilon() );
}

BOOST_FIXTURE_TEST_CASE( monotonic_timestep, clock_fixture )
{
    BOOST_CHECK_NO_THROW( clock.set_timestep(0.001) );
    for (int i = 0; i < 1000000; ++i) {
        clock.advance();
    }
    BOOST_CHECK_EQUAL( clock_const.step(), 1000000u );
    BOOST_CHECK_CLOSE_FRACTION( clock_const.time(), 1000., epsilon() );
}
