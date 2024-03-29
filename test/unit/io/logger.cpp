/*
 * Copyright © 2011  Felix Höfling and Peter Colberg
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

#define BOOST_TEST_MODULE logger
#include <boost/test/unit_test.hpp>

#include <boost/foreach.hpp>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <stdio.h>

#include <halmd/io/logger.hpp>
#include <test/tools/ctest.hpp>

using namespace boost;
using namespace halmd;
using namespace std;

//! count lines of a text file
int count_lines(istream& is)
{
    string line;
    int count = 0;

    while (!is.eof()) {
        if (getline(is, line))
            count++;
        else
            break;
    }
    return count;
}

//! capture console output
class console
{
public:
    //! begin capture to string stream
    console() : buf_(clog.rdbuf(stream_.rdbuf())) {}

    //! end capture to string stream
    ~console()
    {
        clog.rdbuf(buf_);
    }

    //! get string stream
    stringstream& stream()
    {
        return stream_;
    }

private:
    //! string stream
    stringstream stream_;
    //! original console stream buffer
    streambuf* buf_;
};

//! check that logger instance is a singleton
BOOST_AUTO_TEST_CASE( instance )
{
    logging& instance1 = logging::get();
    logging& instance2 = logging::get();
    BOOST_CHECK_EQUAL( &instance1, &instance2 );
}

//! check default log level
BOOST_FIXTURE_TEST_CASE( default_log_level, console )
{
#ifdef NDEBUG
    LOG("message");
    LOG_DEBUG("debug");
    BOOST_CHECK_EQUAL( count_lines(console::stream()), 1 );
#else
    LOG_DEBUG("debug");
    LOG_TRACE("trace");
    BOOST_CHECK_EQUAL( count_lines(console::stream()), 1 );
#endif
}

//! check for duplicate output
BOOST_AUTO_TEST_CASE( duplicate_output )
{
    logging& log = logging::get();
    log.close_console();
    {
        console cons;
        LOG("message");
        BOOST_CHECK_EQUAL( count_lines(cons.stream()), 0 );
    }
    log.open_console(logging::message);
    {
        console cons;
        LOG("message");
        BOOST_CHECK_EQUAL( count_lines(cons.stream()), 1 );
    }
    log.open_console(logging::message);
    {
        console cons;
        LOG("message");
        BOOST_CHECK_EQUAL( count_lines(cons.stream()), 1 );
    }
}

//! test log levels
BOOST_AUTO_TEST_CASE( log_levels )
{
    string const logfile("test_unit_io_logger.log");

    vector<logging::severity_level> log_levels = {
        logging::error
      , logging::warning
      , logging::message
      , logging::info
      , logging::debug
      , logging::trace
    };

    // sweep all verbosities
    BOOST_FOREACH(logging::severity_level log_level, log_levels) {
        logging::get().open_console(log_level);
        logging::get().open_file(logfile, log_level);

        // capture console output
        console cons;

        LOG("message");
        LOG_ERROR("error");
        LOG_WARNING("warning");
        LOG_INFO("info");
        LOG_DEBUG("debugging");
        LOG_TRACE("trace");

        fstream file(logfile.c_str());
#ifdef NDEBUG
        logging::severity_level max_level = min(log_level, logging::info);
        BOOST_CHECK_EQUAL( count_lines(file), max_level + 1 );
        BOOST_CHECK_EQUAL( count_lines(cons.stream()), max_level + 1 );
#else
        BOOST_CHECK_EQUAL( count_lines(file), log_level + 1 );
        BOOST_CHECK_EQUAL( count_lines(cons.stream()), log_level + 1 );
#endif

        logging::get().close_console();
        logging::get().close_file();
    }

    remove(logfile.c_str());
}
