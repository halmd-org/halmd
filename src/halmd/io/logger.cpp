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

#include <boost/log/attributes/clock.hpp>
#include <boost/log/filters/attr.hpp>
#include <boost/log/formatters/attr.hpp>
#include <boost/log/formatters/date_time.hpp>
#include <boost/log/formatters/format.hpp>
#include <boost/log/formatters/message.hpp>
#include <boost/log/utility/empty_deleter.hpp>

#include <halmd/io/logger.hpp>

using namespace boost;
using namespace boost::log;
using namespace std;

#define TIMESTAMP_FORMAT "%d-%m-%Y %H:%M:%S.%f"

namespace halmd
{
namespace io
{

/**
 * Assemble module options
 */
void logging::options(po::options_description& desc)
{
    desc.add_options()
        ("verbose,v", po::accum_value<int>()->default_value(warning),
         "increase verbosity")
        ;
}

sources::severity_logger<logging::severity_level> logging::logger;

static inline ostream& operator<<(ostream& os, logging::severity_level level)
{
    switch (level)
    {
      case logging::trace:
        os << "[TRACE] "; break;
      case logging::debug:
        os << "[DEBUG] "; break;
      case logging::warning:
        os << "[WARNING] "; break;
      case logging::error:
        os << "[ERROR] "; break;
      case logging::fatal:
        os << "[FATAL] "; break;
      default:
        break;
    }
    return os;
}

logging::logging(po::options const& vm)
{
    {
        shared_ptr<console_backend> backend(make_shared<console_backend>());
        backend->add_stream(
            shared_ptr<ostream>(&clog, empty_deleter())
        );
        backend->set_formatter(
            formatters::format("[%1%] %2%%3%")
                % formatters::date_time("TimeStamp", keywords::format = TIMESTAMP_FORMAT)
                % formatters::attr<severity_level>("Severity")
                % formatters::message()
        );
        backend->auto_flush(true);

        console_ = make_shared<console_sink>(backend);
        console_->set_filter(
            filters::attr<severity_level>("Severity") <= vm["verbose"].as<int>()
        );
        core::get()->add_sink(console_);
    }

    {
        shared_ptr<file_backend> backend(
            make_shared<file_backend>(
                keywords::file_name = vm["output"].as<string>() + ".log"
            )
        );
        backend->set_formatter(
            formatters::format("[%1%] %2%%3%")
                % formatters::date_time("TimeStamp", keywords::format = TIMESTAMP_FORMAT)
                % formatters::attr<severity_level>("Severity")
                % formatters::message()
        );
        backend->auto_flush(true);

        file_ = make_shared<file_sink>(backend);
        file_->set_filter(
            filters::attr<severity_level>("Severity") <= max(
                vm["verbose"].as<int>()
              , static_cast<int>(info)
            )
        );
        core::get()->add_sink(file_);
    }

    core::get()->add_global_attribute(
        "TimeStamp"
      , make_shared<attributes::local_clock>()
    );
}

logging::~logging()
{
    core::get()->remove_sink(console_);
    core::get()->remove_sink(file_);
}

} // namespace io

} // namespace halmd
