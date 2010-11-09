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

#define BOOST_TEST_MODULE test_io_logger
#include <boost/test/unit_test.hpp>

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

struct timer_map {
    // define and register profiling timers
    HALMD_PROFILING_TAG( timer1, "first timer" );
    HALMD_PROFILING_TAG( timer2, "second timer" );

    boost::fusion::map<
        fusion::pair<timer1, accumulator<double> >
      , fusion::pair<timer2, accumulator<double> >
    > map;
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

    // resolve module dependencies
    using namespace halmd::io::profiling;
    shared_ptr<utility::profiler> profiler;

    // repeat three times
    for (unsigned n=0; n < 3; n++) {
        BOOST_TEST_MESSAGE("Pass #" << n+1);

        // construct modules
        vector<shared_ptr<writer> > writers;
        writers.push_back(make_shared<writers::log>());
        writers.push_back(make_shared<writers::hdf5>(file_name));
        profiler = make_shared<utility::profiler>(writers);

        // register profiling timers
        timer_map timers;
        profiler->register_map(timers.map);

        // accumulate some values
        for (float x=0; x < 1; x += 0.1) {
            fusion::at_key<timer_map::timer1>(timers.map)(x);
            fusion::at_key<timer_map::timer2>(timers.map)(n * x);
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
