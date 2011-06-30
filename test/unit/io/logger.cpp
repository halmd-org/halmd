/*
 * Copyright © 2011  Felix Höfling and Peter Colberg
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

#define BOOST_TEST_MODULE logger
#include <boost/test/unit_test.hpp>

#include <boost/shared_ptr.hpp>
#include <fstream>
#include <iostream>
#include <map>
#include <stdio.h>

#include <halmd/io/logger.hpp>

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
    logger& instance1 = logger::instance();
    logger& instance2 = logger::instance();
    BOOST_CHECK_EQUAL( &instance1, &instance2 );
}

//! check default log level
BOOST_FIXTURE_TEST_CASE( default_log_level, console )
{
#ifdef NDEBUG
    LOG("info");
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
    logger& log = logger::instance();
    log.close_console();
    {
        console cons;
        LOG("info");
        BOOST_CHECK_EQUAL( count_lines(cons.stream()), 0 );
    }
    log.open_console(logger::info);
    {
        console cons;
        LOG("info");
        BOOST_CHECK_EQUAL( count_lines(cons.stream()), 1 );
    }
    log.open_console(logger::info);
    {
        console cons;
        LOG("info");
        BOOST_CHECK_EQUAL( count_lines(cons.stream()), 1 );
    }
}

//! test log levels
BOOST_AUTO_TEST_CASE( log_levels )
{
    string const logfile("test_unit_io_logger.log");

    // sweep all verbosities
    for (int log_level = -1; log_level <= 5; log_level++) {
        logger::severity_level s = static_cast<logger::severity_level>(log_level);
        logger::instance().open_console(s);
        logger::instance().open_file(logfile, s);

        // capture console output
        console cons;

        LOG("info");
        LOG_FATAL("fatal");
        LOG_ERROR("error");
        LOG_WARNING("warning");
        LOG_DEBUG("debugging");
        LOG_TRACE("trace");

        fstream file(logfile.c_str());
#ifdef NDEBUG
        int max_level = min(log_level, static_cast<int>(logger::info));
        BOOST_CHECK_EQUAL( count_lines(file), max_level + 1 );
        BOOST_CHECK_EQUAL( count_lines(cons.stream()), max_level + 1 );
#else
        BOOST_CHECK_EQUAL( count_lines(file), log_level + 1 );
        BOOST_CHECK_EQUAL( count_lines(cons.stream()), log_level + 1 );
#endif

        logger::instance().close_console();
        logger::instance().close_file();
    }

    remove(logfile.c_str());
}
