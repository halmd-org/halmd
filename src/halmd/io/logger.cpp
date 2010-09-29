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
#include <halmd/utility/luabind.hpp>

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
        ("verbose,v", accum_value<int>()->default_value(warning),
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

logging::logging(po::variables_map const& vm)
{
    core::get()->add_global_attribute(
        "TimeStamp"
      , make_shared<attributes::local_clock>()
    );
    log_to_console(
        static_cast<severity_level>(vm["verbose"].as<int>())
    );
    log_to_file(
        static_cast<severity_level>(
            max(vm["verbose"].as<int>(), static_cast<int>(info))
        )
      , vm["output"].as<string>() + ".log"
    );
}

/**
 * enable logging to console
 *
 * @param level logging severity level
 *
 * FIXME repeated calls of this function if public
 */
void logging::log_to_console(severity_level level)
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
        filters::attr<severity_level>("Severity") <= level
    );
    core::get()->add_sink(console_);
}

/**
 * enable logging to file
 *
 * @param level logging severity level
 *
 * FIXME repeated calls of this function if public
 */
void logging::log_to_file(severity_level level, string file_name)
{
    shared_ptr<file_backend> backend(
        make_shared<file_backend>(
            keywords::file_name = file_name
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
        filters::attr<severity_level>("Severity") <= level
    );
    core::get()->add_sink(file_);
}

/**
 * remove sinks from logging core singleton
 */
logging::~logging()
{
    core::get()->remove_sink(console_);
    core::get()->remove_sink(file_);
}

template <enum logging::severity_level Level>
static void log_wrapper(char const* message)
{
    BOOST_LOG_SEV(logging::logger, Level) << message;
}

static __attribute__((constructor)) void register_lua()
{
    using namespace luabind;
    lua_registry::get()->push_back
    ((
        namespace_("halmd_wrapper")
        [
            namespace_("io")
            [
                class_<logging, shared_ptr<logging> >("logging")
                    .scope
                    [
                        def("options", &logging::options)
                    ]
            ]
        ]
      , namespace_("log")
        [
            def("fatal", &log_wrapper<logging::fatal>)
          , def("error", &log_wrapper<logging::error>)
          , def("warning", &log_wrapper<logging::warning>)
          , def("info", &log_wrapper<logging::info>)
          , def("debug", &log_wrapper<logging::debug>)
          , def("trace", &log_wrapper<logging::trace>)
        ]
    ));
}

} // namespace io

} // namespace halmd
