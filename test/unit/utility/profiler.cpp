/*
 * Copyright © 2010  Felix Höfling
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

#define BOOST_TEST_MODULE profiler
#include <boost/test/unit_test.hpp>

#include <memory>

#include <halmd/utility/profiler.hpp>
#include <test/tools/ctest.hpp>

using namespace boost;
using namespace halmd;
using namespace std;

/**
 * profiling timers
 */
struct timer_map
{
    typedef utility::profiler::accumulator_type accumulator_type;
    std::shared_ptr<accumulator_type> timer1;
    std::shared_ptr<accumulator_type> timer2;
};

//
// test profiler module
//

BOOST_AUTO_TEST_CASE( test_profiler )
{
    std::shared_ptr<utility::profiler> profiler;
    timer_map timers;
    timers.timer1 = std::make_shared<timer_map::accumulator_type>();
    timers.timer2 = std::make_shared<timer_map::accumulator_type>();

    // repeat three times
    for (unsigned n=0; n < 3; n++) {
        BOOST_TEST_MESSAGE("Pass #" << n+1);

        // construct modules
        profiler = std::make_shared<utility::profiler>();

        // register profiling timers
        profiler->on_profile(timers.timer1, "first timer");
        profiler->on_profile(timers.timer2, "second timer");

        // accumulate some values
        for (float x=0; x < 1; x += 0.1) {
            (*timers.timer1)(x);
            (*timers.timer2)(n * x);
        }

        profiler->profile();
    }

    // FIXME add some tests here (e.g. line counting of log file)
}
