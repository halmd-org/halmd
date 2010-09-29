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

#undef NDEBUG
#include <halmd/io/logger.hpp>
#define NDEBUG
#include <halmd/options.hpp>


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
// test logger module
//

BOOST_AUTO_TEST_CASE( test_logger )
{
    typedef boost::program_options::variable_value variable_value;

    string const logfile("halmd_test.log");

    // sweep all verbosities
    for (int log_level = -1; log_level <= 5; log_level++) {

        shared_ptr<logger> log(new logger);
        log->log_to_console((logger::severity_level)log_level);
        log->log_to_file((logger::severity_level)log_level, logfile);

        LOG("info");
        LOG_FATAL("fatal");
        LOG_ERROR("error");
        LOG_WARNING("warning");
        LOG_DEBUG("debugging");
        LOG_TRACE("trace");
        BOOST_CHECK_EQUAL((int)count_lines(logfile.c_str()), log_level + 1);

        // drop logger and close file
        log.reset();
    }

    remove(logfile.c_str());
}

