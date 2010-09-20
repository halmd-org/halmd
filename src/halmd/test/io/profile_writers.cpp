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
#include <halmd/io/profile/writers/log.hpp>
#include <halmd/io/profile/writers/hdf5.hpp>
#include <halmd/utility/options.hpp>
#include <halmd/utility/profiler.hpp>

#include <halmd/utility/module.hpp>
#include <halmd/utility/modules/factory.hpp>
#include <halmd/utility/modules/policy.hpp>
#include <halmd/utility/modules/resolver.hpp>

using namespace boost;
using namespace halmd;
using namespace std;

struct timer_map {
    // define and register profile timers
    HALMD_PROFILE_TAG( timer1, "first timer" );
    HALMD_PROFILE_TAG( timer2, "second timer" );

    boost::fusion::map<
        fusion::pair<timer1, accumulator<double> >
      , fusion::pair<timer2, accumulator<double> >
    > map;
};

//
// test profile writer modules
//

BOOST_AUTO_TEST_CASE( test_profile_writers )
{
    typedef boost::program_options::variable_value variable_value;

    halmd::po::options vm;
    // override const operator[] in variables_map
    map<string, variable_value>& vm_(vm);
    vm_["verbose"] = variable_value(0, true);
    vm_["output"] = variable_value(string("halmd_test"), true);

    // enable logging to console
    shared_ptr<logging> logger(new logging(vm));

    // resolve module dependencies
    using namespace halmd::io::profile;
    po::unparsed_options unparsed;
    modules::resolver resolver(modules::registry::graph());
    resolver.resolve<writers::log>(vm, unparsed);
    resolver.resolve<writers::hdf5>(vm, unparsed);
    resolver.resolve<utility::profiler>(vm, unparsed);
    modules::policy policy(resolver.graph());
    modules::factory factory(policy.graph());

    shared_ptr<writers::log> log;
    shared_ptr<writers::hdf5> hdf5;
    shared_ptr<utility::profiler> profiler;

    // repeat three times
    for (unsigned n=0; n < 3; n++) {
        BOOST_TEST_MESSAGE("Pass #" << n+1);

        // construct modules
        log = modules::fetch<writers::log>(factory, vm);
        hdf5 = modules::fetch<writers::hdf5>(factory, vm);
        profiler = modules::fetch<utility::profiler>(factory, vm);

        // register profile timers
        timer_map timers;
        profiler->register_map(timers.map);

        // accumulate some values
        for (float x=0; x < 1; x += 0.1) {
            fusion::at_key<timer_map::timer1>(timers.map)(x);
            fusion::at_key<timer_map::timer2>(timers.map)(n * x);
        }

        // write results
        hdf5->write();
        log->write();

        // destroy some modules
        log.reset();
        if (n < 1)
            hdf5.reset();
        // if the profiler is not destroyed as well, the test fails:
        // HDF5-DIAG: Error detected in HDF5 (1.8.1) thread 0:
        // #000: H5D.c line 171 in H5Dcreate2(): unable to create dataset
        if (n < 2)
            profiler.reset();
    }

    // FIXME add some tests here (e.g. line counting of log file)

    // remove files
    remove((vm["output"].as<string>() + ".prf").c_str());
    remove((vm["output"].as<string>() + ".log").c_str());
}

