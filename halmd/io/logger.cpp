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

// increase compiler compatibility, e.g. with Clang 2.8
#define BOOST_LOG_NO_UNSPECIFIED_BOOL
#include <boost/log/attributes/clock.hpp>
#include <boost/log/filters/attr.hpp>
#include <boost/log/formatters/attr.hpp>
#include <boost/log/formatters/date_time.hpp>
#include <boost/log/formatters/format.hpp>
#include <boost/log/formatters/message.hpp>
#include <boost/log/utility/empty_deleter.hpp>
#include <boost/version.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace boost::log;
using namespace std;

#define TIMESTAMP_FORMAT "%d-%m-%Y %H:%M:%S.%f"

namespace halmd
{

static inline ostream& operator<<(ostream& os, logger::severity_level level)
{
    switch (level)
    {
      case logger::trace:
        os << "[TRACE] "; break;
      case logger::debug:
        os << "[DEBUG] "; break;
      case logger::warning:
        os << "[WARNING] "; break;
      case logger::error:
        os << "[ERROR] "; break;
      case logger::fatal:
        os << "[FATAL] "; break;
      default:
        break;
    }
    return os;
}

logger::logger()
{
    core::get()->add_global_attribute(
        "TimeStamp"
#ifdef BOOST_LOG_ATTRIBUTE_HPP_INCLUDED_ // Boost.Log < r479 (SVN)
      , make_shared<attributes::local_clock>()
#else
      , attributes::local_clock()
#endif
    );
}

/**
 * enable logger to console
 *
 * @param level logger severity level
 *
 * FIXME repeated calls of this function if public
 */
void logger::log_to_console(severity_level level)
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

    core::get()->remove_sink(console_);
    console_ = make_shared<console_sink>(backend);
    console_->set_filter(
        filters::attr<severity_level>("Severity") <= level
    );
    core::get()->add_sink(console_);
}

/**
 * enable logger to file
 *
 * @param level logger severity level
 *
 * FIXME repeated calls of this function if public
 */
void logger::log_to_file(severity_level level, string file_name)
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

    core::get()->remove_sink(file_);
    file_ = make_shared<file_sink>(backend);
    file_->set_filter(
        filters::attr<severity_level>("Severity") <= level
    );
    core::get()->add_sink(file_);
}

/**
 * remove sinks from logger core singleton
 */
logger::~logger()
{
    core::get()->remove_sink(console_);
    core::get()->remove_sink(file_);
}

sources::severity_logger<logger::severity_level> logger::logger_;

template <enum logger::severity_level level>
static void log_wrapper(char const* message)
{
    BOOST_LOG_SEV(logger::get(), level) << message;
}

void logger::luaopen(lua_State* L)
{
    using namespace luabind;
    module(L)
    [
        namespace_("libhalmd")
        [
            namespace_("io")
            [
                namespace_("logger")
                [
                    def("fatal", &log_wrapper<logger::fatal>)
                  , def("error", &log_wrapper<logger::error>)
                  , def("warning", &log_wrapper<logger::warning>)
                  , def("info", &log_wrapper<logger::info>)
#ifndef NDEBUG
                  , def("debug", &log_wrapper<logger::debug>)
                  , def("trace", &log_wrapper<logger::trace>)
#endif
                ]
            ]
        ]
    ];
}

HALMD_LUA_API int luaopen_libhalmd_logger(lua_State* L)
{
    logger::luaopen(L);
    return 0;
}

} // namespace halmd
