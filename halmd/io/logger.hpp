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

#ifndef HALMD_IO_LOGGER_HPP
#define HALMD_IO_LOGGER_HPP

// increase compiler compatibility, e.g. with Clang 2.8
#define BOOST_LOG_NO_UNSPECIFIED_BOOL
#include <boost/log/sinks/sync_frontend.hpp>
#include <boost/log/sinks/text_file_backend.hpp>
#include <boost/log/sinks/text_ostream_backend.hpp>
#include <boost/log/sources/record_ostream.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/noncopyable.hpp>
#include <lua.hpp>

namespace halmd {

//! Logger
class logger
  : boost::noncopyable
{
public:
    //! logging severity levels
    enum severity_level
    {
        fatal,
        error,
        warning,
        info,
        debug,
        trace,
    };
    //! logger source type
    typedef boost::log::sources::severity_logger<severity_level> logger_type;

    //! open log to console
    void open_console(severity_level level);
    //! close log to console
    void close_console();
    //! open log to file
    void open_file(std::string file_name, severity_level level);
    //! close log to file
    void close_file();
    //! Lua bindings
    static void luaopen(lua_State* L);

    //! get logger source
    logger_type& get()
    {
       return logger_;
    }

    //! get logger singleton instance
    static logger& instance()
    {
        static logger logger_;
        return logger_;
    }

private:
    typedef boost::log::sinks::text_ostream_backend console_backend_type;
    typedef boost::log::sinks::synchronous_sink<console_backend_type> console_sink_type;
    typedef boost::log::sinks::text_file_backend file_backend_type;
    typedef boost::log::sinks::synchronous_sink<file_backend_type> file_sink_type;

    //! initialize logger
    logger();

    //! log timestamp format
    static char const* timestamp()
    {
        return "%d-%m-%Y %H:%M:%S.%f";
    }

    //! logger source
    boost::log::sources::severity_logger<severity_level> logger_;
    //! console sink
    boost::shared_ptr<console_sink_type> console_;
    //! file sink
    boost::shared_ptr<file_sink_type> file_;
};

#define HALMD_LOG(level, format)                        \
{                                                       \
    BOOST_LOG_SEV(                                      \
        halmd::logger::instance().get()                 \
      , level                                           \
    ) << format;                                        \
}

#define HALMD_LOG_ONCE(level, format)                   \
{                                                       \
    static bool logged = false;                         \
    if (!logged) {                                      \
        BOOST_LOG_SEV(                                  \
            halmd::logger::instance().get()             \
          , level                                       \
        ) << format;                                    \
        logged = true;                                  \
    }                                                   \
}

#define LOG_FATAL(format)           HALMD_LOG(halmd::logger::fatal, format)
#define LOG_FATAL_ONCE(format)      HALMD_LOG_ONCE(halmd::logger::fatal, format)
#define LOG_ERROR(format)           HALMD_LOG(halmd::logger::error, format)
#define LOG_ERROR_ONCE(format)      HALMD_LOG_ONCE(halmd::logger::error, format)
#define LOG_WARNING(format)         HALMD_LOG(halmd::logger::warning, format)
#define LOG_WARNING_ONCE(format)    HALMD_LOG_ONCE(halmd::logger::warning, format)
#define LOG(format)                 HALMD_LOG(halmd::logger::info, format)
#define LOG_ONCE(format)            HALMD_LOG_ONCE(halmd::logger::info, format)
#ifndef NDEBUG
# define LOG_DEBUG(format)          HALMD_LOG(halmd::logger::debug, format)
# define LOG_DEBUG_ONCE(format)     HALMD_LOG_ONCE(halmd::logger::debug, format)
# define LOG_TRACE(format)          HALMD_LOG(halmd::logger::trace, format)
# define LOG_TRACE_ONCE(format)     HALMD_LOG_ONCE(halmd::logger::trace, format)
#else
# define LOG_DEBUG(format)
# define LOG_DEBUG_ONCE(format)
# define LOG_TRACE(format)
# define LOG_TRACE_ONCE(format)
#endif

} // namespace halmd

#endif /* ! HALMD_IO_LOGGER_HPP */
