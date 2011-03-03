/*
 * Copyright © 2010  Felix Höfling
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

#define BOOST_TEST_MODULE test_profiling_writers
#include <boost/test/unit_test.hpp>

#include <boost/assign.hpp>
#include <boost/shared_ptr.hpp>
#include <fstream>
#include <map>
#include <stdio.h>

#include <halmd/io/logger.hpp>
#include <halmd/io/profiling/writers/log.hpp>
#include <halmd/io/profiling/writers/hdf5.hpp>
#include <halmd/utility/profiler.hpp>

using namespace boost;
using namespace halmd;
using namespace std;

/**
 * profiling timers
 */
struct timer_map
{
    typedef utility::profiler::accumulator_type accumulator_type;
    accumulator_type timer1;
    accumulator_type timer2;
};

//
// test profiling writer modules
//

BOOST_AUTO_TEST_CASE( test_profiling_writers )
{
    string const file_name("test_io_logger.prf");

    // enable logging to console
    static logger log;
    log.log_to_console(logger::debug);

    using namespace halmd::io::profiling;
    using namespace boost::assign;
    shared_ptr<utility::profiler> profiler;

    // repeat three times
    for (unsigned n=0; n < 3; n++) {
        BOOST_TEST_MESSAGE("Pass #" << n+1);

        // construct modules
        vector<shared_ptr<writer> > writers;
        writers.push_back(make_shared<writers::log>());
        writers.push_back(make_shared<writers::hdf5>(file_name));
        vector<string> tag = list_of("test")("timer_map");
        profiler = make_shared<utility::profiler>(writers, tag);

        // register profiling timers
        timer_map timers;
        profiler->register_runtime(timers.timer1, "first timer");
        profiler->register_runtime(timers.timer2, "second timer");

        // accumulate some values
        for (float x=0; x < 1; x += 0.1) {
            timers.timer1(x);
            timers.timer2(n * x);
        }

        // write results
        for_each(writers.begin(), writers.end(), bind(&writer::write, _1));

        // destroy some modules
        writers.pop_back();
        if (n < 1)
            writers.pop_back();
        // if the profiler is not destroyed as well, the test fails:
        // HDF5-DIAG: Error detected in HDF5 (1.8.1) thread 0:
        // #000: H5D.c line 171 in H5Dcreate2(): unable to create dataset
        if (n < 2)
            profiler.reset();
    }

    // FIXME add some tests here (e.g. line counting of log file)

    // remove files
    remove(file_name.c_str());
}
