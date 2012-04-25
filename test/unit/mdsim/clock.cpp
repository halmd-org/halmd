/*
 * Copyright © 2011  Peter Colberg
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

#define BOOST_TEST_MODULE clock
#include <boost/test/unit_test.hpp>

#include <boost/make_shared.hpp>
#include <limits>

#include <halmd/mdsim/clock.hpp>
#include <test/tools/ctest.hpp>

using namespace boost;
using namespace std;

typedef halmd::mdsim::clock clock_type;
typedef clock_type::step_type step_type;
typedef clock_type::time_type time_type;

time_type const eps_time = numeric_limits<time_type>::epsilon();

BOOST_AUTO_TEST_CASE( init )
{
    shared_ptr<clock_type> clock = make_shared<clock_type>(0.001);
    BOOST_CHECK_EQUAL( clock->step(), 0u );
    BOOST_CHECK_CLOSE_FRACTION( clock->time(), 0., eps_time );
    BOOST_CHECK_CLOSE_FRACTION( clock->timestep(), 0.001, eps_time );
}

BOOST_AUTO_TEST_CASE( monotonic_timestep )
{
    shared_ptr<clock_type> clock = make_shared<clock_type>(0.001);
    for (int i = 0; i < 1000000; ++i) {
        clock->advance();
    }
    BOOST_CHECK_EQUAL( clock->step(), 1000000u );
    BOOST_CHECK_CLOSE_FRACTION( clock->time(), 1000., eps_time );
}
