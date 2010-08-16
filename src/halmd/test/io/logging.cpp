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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_io_logger
#include <boost/test/unit_test.hpp>

#include <boost/shared_ptr.hpp>
#include <fstream>
#include <map>
#include <stdio.h>

#undef NDEBUG
#include <halmd/io/logger.hpp>
#define NDEBUG
#include <halmd/utility/options.hpp>


using namespace boost;
using namespace halmd;
using namespace std;

// count lines of a text file
unsigned count_lines(char const* filename)
{
    ifstream f(filename);
    string line;
    unsigned count = 0;

    while (!f.eof()) {
        if (getline(f, line))
            count++;
        else
            break;
    }
    return count;
}

//
// test logging module
//

BOOST_AUTO_TEST_CASE( test_logging )
{
    typedef boost::program_options::variable_value variable_value;

    halmd::po::options vm;
    // override const operator[] in variables_map
    map<string, variable_value>& vm_(vm);
    vm_["verbose"] = variable_value(0, true);
    vm_["output"] = variable_value(string("test_logging"), true);

    string logfile = vm["output"].as<string>() + ".log";

    // sweep all verbosities
    for (int log_level = -1; log_level <= 5; log_level++) {

        vm_["verbose"] = variable_value(log_level, false);
        shared_ptr<logging> logger(new logging(vm));

        LOG("info");
        LOG_FATAL("fatal");
        LOG_ERROR("error");
        LOG_WARNING("warning");
        LOG_DEBUG("debugging");
        LOG_TRACE("trace");
        BOOST_CHECK_EQUAL((int)count_lines(logfile.c_str()), max(4, log_level + 1));

        // drop logger and close file
        logger.reset();
    }

    remove(logfile.c_str());
}

