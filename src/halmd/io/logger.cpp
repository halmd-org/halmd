/*
 * Copyright © 2008-2010  Peter Colberg
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

#include <boost/log/attributes/timer.hpp>
#include <boost/log/filters/attr.hpp>
#include <boost/log/formatters/attr.hpp>
#include <boost/log/formatters/date_time.hpp>
#include <boost/log/formatters/format.hpp>
#include <boost/log/formatters/message.hpp>
#include <boost/log/utility/init/common_attributes.hpp>
#include <boost/log/utility/init/to_console.hpp>
#include <boost/log/utility/init/to_file.hpp>

#include <halmd/io/logger.hpp>

using namespace boost::log;
using namespace boost;
using namespace std;

#define TIMESTAMP_FORMAT "%d-%m-%Y %H:%M:%S.%f"

namespace halmd
{

sources::severity_logger<severity_level> logger_;

static inline ostream& operator<<(ostream& os, severity_level lvl)
{
    switch (lvl)
    {
      case trace:
        os << "[TRACE] "; break;
      case debug:
        os << "[DEBUG] "; break;
      case warning:
        os << "[WARNING] "; break;
      case error:
        os << "[ERROR] "; break;
      case fatal:
        os << "[FATAL] "; break;
      default:
        break;
    }
    return os;
}

namespace io { namespace logger
{

/**
 * initialize logging
 */
void init(options const& vm)
{
    severity_level lvl_cons, lvl_file;
    switch (vm["verbose"].as<int>())
    {
      case 0:
        lvl_cons = warning; lvl_file = info; break;
      case 1:
        lvl_cons = lvl_file = info; break;
      case 2:
        lvl_cons = lvl_file = debug; break;
      default:
        lvl_cons = lvl_file = trace; break;
    }

    init_log_to_file
    (
        vm["output"].as<string>(),
        keywords::auto_flush = true,
        keywords::filter = filters::attr<severity_level>("Severity") >= lvl_file,
        keywords::format = formatters::format("[%1%] %2%%3%")
            % formatters::date_time("TimeStamp", keywords::format = TIMESTAMP_FORMAT)
            % formatters::attr<severity_level>("Severity")
            % formatters::message()
    );

    init_log_to_console
    (
        std::clog,
        keywords::auto_flush = true,
        keywords::filter = filters::attr<severity_level>("Severity") >= lvl_cons,
        keywords::format = formatters::format("[%1%] %2%%3%")
            % formatters::date_time("TimeStamp", keywords::format = TIMESTAMP_FORMAT)
            % formatters::attr<severity_level>("Severity")
            % formatters::message()
    );

    // add timestamp and record counter
    add_common_attributes();
}

}} // namespace io::logger

} // namespace halmd
