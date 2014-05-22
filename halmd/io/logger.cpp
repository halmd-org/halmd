/*
 * Copyright © 2018 Felix Höfling
 * Copyright © 2014 Nicolas Höft
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
#include <boost/log/attributes/constant.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/support/date_time.hpp>
#include <boost/utility/empty_deleter.hpp>
#include <boost/version.hpp>

#include <halmd/io/logger.hpp>
#include <halmd/utility/lua/lua.hpp>

using namespace boost::log;

namespace halmd {

BOOST_LOG_ATTRIBUTE_KEYWORD(severity, "Severity", logging::severity_level)
BOOST_LOG_ATTRIBUTE_KEYWORD(module_attr, "Module", std::string)
BOOST_LOG_ATTRIBUTE_KEYWORD(tag_attr, "Tag", std::string)

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
    boost::shared_ptr<console_backend_type> backend(boost::make_shared<console_backend_type>());
    backend->add_stream(
        boost::shared_ptr<std::ostream>(&std::clog, boost::empty_deleter())
    );
    backend->auto_flush(true);

    core::get()->remove_sink(console_);
    console_ = boost::make_shared<console_sink_type>(backend);
    console_->set_filter(
#ifdef NDEBUG
        severity <= std::min(level, info)
#else
        severity <= level
#endif
    );
    set_formatter(console_);
    core::get()->add_sink(console_);
}

void logging::close_console()
{
    core::get()->remove_sink(console_);
    console_.reset();
}

void logging::open_file(std::string file_name, severity_level level)
{
    boost::shared_ptr<file_backend_type> backend(
        boost::make_shared<file_backend_type>(
            keywords::file_name = file_name
        )
    );
    backend->auto_flush(true);

    core::get()->remove_sink(file_);
    file_ = boost::make_shared<file_sink_type>(backend);
    file_->set_filter(
#ifdef NDEBUG
        severity <= std::min(level, info)
#else
        severity <= level
#endif
    );
    set_formatter(file_);
    core::get()->add_sink(file_);
}

void logging::close_file()
{
    core::get()->remove_sink(file_);
    file_.reset();
}

static inline std::ostream& operator<<(std::ostream& os, logging::severity_level level)
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
void logging::set_formatter(boost::shared_ptr<backend_type> backend)
{
    backend->set_formatter(expressions::stream
        << expressions::format_date_time<boost::posix_time::ptime>("TimeStamp", "[%d-%m-%Y %H:%M:%S.%f]")
        << expressions::if_(severity != logging::info)
           [
               expressions::stream << " [" << severity << "]"
           ]
        << expressions::if_(expressions::has_attr(module_attr))
           [
               expressions::stream
                   << " " << module_attr
                   << expressions::if_(expressions::has_attr(tag_attr))
                      [
                          expressions::stream << "[" << tag_attr << "]"
                      ]
                   << ":"
           ]
        << " " << expressions::smessage
    );
}

template <logging::severity_level level>
static void wrap_log(char const* message)
{
    HALMD_LOG(level, message);
}

template <typename T>
static void wrap_add_attribute(logger& logger_, std::string const& attr, T const& value)
{
    logger_.add_attribute(attr, boost::log::attributes::constant<T>(value));
}

void logging::luaopen(lua_State* L)
{
    using namespace luabind;
    module(L, "libhalmd")
    [
        namespace_("io")
        [
            class_<logger, boost::shared_ptr<logger> >("logger")
                .def(constructor<>())
                .def("add_attribute", &wrap_add_attribute<std::string>)
                .def("add_attribute", &wrap_add_attribute<int>)

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
boost::shared_ptr<logger> const logger_ = boost::make_shared<logger>();

HALMD_LUA_API int luaopen_libhalmd_io_logger(lua_State* L)
{
    logging::luaopen(L);
    return 0;
}

} // namespace halmd
