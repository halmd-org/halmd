/*
 * Copyright © 2008-2011  Peter Colberg
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
#include <boost/log/formatters/if.hpp>
#include <boost/log/formatters/message.hpp>
#include <boost/log/formatters/stream.hpp>
#include <boost/log/utility/empty_deleter.hpp>
#include <boost/version.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost;
using namespace boost::log;
using namespace std;

namespace halmd {

logging::logging()
{
    core::get()->add_global_attribute("TimeStamp", attributes::local_clock());
#ifdef NDEBUG
    logging::open_console(info);
#else
    logging::open_console(debug);
#endif
}

void logging::open_console(severity_level level)
{
    shared_ptr<console_backend_type> backend(make_shared<console_backend_type>());
    backend->add_stream(
        shared_ptr<ostream>(&clog, empty_deleter())
    );
    set_formatter(backend);
    backend->auto_flush(true);

    core::get()->remove_sink(console_);
    console_ = make_shared<console_sink_type>(backend);
    console_->set_filter(
        filters::attr<severity_level>("Severity") <= level
    );
    core::get()->add_sink(console_);
}

void logging::close_console()
{
    core::get()->remove_sink(console_);
    console_.reset();
}

void logging::open_file(string file_name, severity_level level)
{
    shared_ptr<file_backend_type> backend(
        make_shared<file_backend_type>(
            keywords::file_name = file_name
        )
    );
    set_formatter(backend);
    backend->auto_flush(true);

    core::get()->remove_sink(file_);
    file_ = make_shared<file_sink_type>(backend);
    file_->set_filter(
        filters::attr<severity_level>("Severity") <= level
    );
    core::get()->add_sink(file_);
}

void logging::close_file()
{
    core::get()->remove_sink(file_);
    file_.reset();
}

static inline ostream& operator<<(ostream& os, logging::severity_level level)
{
    switch (level)
    {
      case logging::trace:
        os << "TRACE"; break;
      case logging::debug:
        os << "DEBUG"; break;
      case logging::warning:
        os << "WARNING"; break;
      case logging::error:
        os << "ERROR"; break;
      case logging::fatal:
        os << "FATAL"; break;
      default:
        os << static_cast<int>(level); break;
    }
    return os;
}

template <typename backend_type>
void logging::set_formatter(shared_ptr<backend_type> backend) const
{
    backend->set_formatter(formatters::stream
        << formatters::date_time("TimeStamp", "[%d-%m-%Y %H:%M:%S.%f]")
        << formatters::if_(filters::attr<severity_level>("Severity") != logging::info)
           [
               formatters::stream << " [" << formatters::attr<logging::severity_level>("Severity") << "]"
           ]
        << " " << formatters::message()
    );
}

template <logging::severity_level level>
static void wrap_log(char const* message)
{
    HALMD_LOG(level, message);
}

void logging::luaopen(lua_State* L)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("io")
        [
            class_<logger, shared_ptr<logger> >("logger")

          , namespace_("logging")
            [
                def("fatal", &wrap_log<logging::fatal>)
              , def("error", &wrap_log<logging::error>)
              , def("warning", &wrap_log<logging::warning>)
              , def("info", &wrap_log<logging::info>)
#ifndef NDEBUG
              , def("debug", &wrap_log<logging::debug>)
              , def("trace", &wrap_log<logging::trace>)
#endif
            ]
        ]
    ];
}

/** define logging singleton instance */
logging logging::logging_;
/** define global logger source */
shared_ptr<logger> const logger_ = make_shared<logger>();

HALMD_LUA_API int luaopen_libhalmd_io_logger(lua_State* L)
{
    logging::luaopen(L);
    return 0;
}

} // namespace halmd
