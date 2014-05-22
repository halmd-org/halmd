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
#include <boost/shared_ptr.hpp>
#include <lua.hpp>

namespace halmd {

/**
 * Logging to console and file
 *
 * This class implements a logging interface with severity levels.
 * More accurately it adds a console and a file sink to the log core,
 * and supports discarding messages below a desired severity level.
 *
 * We use the Boost.Log library, which differentiates between log
 * sources and log sinks. The logging class handles only log sinks,
 * and therefore is needed only to setup logging console and file.
 *
 * http://boost-log.sourceforge.net/
 */
class logging
  : boost::noncopyable
{
public:
    /**
     * Log severity levels
     *
     * Note that lower severity levels follow *after* higher severity
     * levels, which makes it more convenient to parse the severity
     * level from an integer option value.
     */
    enum severity_level
    {
        fatal,
        error,
        warning,
        info,
        debug,
        trace
    };

    /**
     * open log to console with given severity level
     *
     * This method may be called repeatedly to set a different severity level.
     */
    void open_console(severity_level level);
    /** close log to console */
    void close_console();
    /**
     * open log to file with given severity level
     *
     * This method may be called repeatedly to set a different severity level.
     */
    void open_file(std::string file_name, severity_level level);
    /** close log to file */
    void close_file();
    /** Lua bindings */
    static void luaopen(lua_State* L);

    /**
     * get logger singleton instance
     *
     * To ensure that logging is enabled by default, which is especially
     * convenient in the unit tests, we make logging a singleton instance.
     */
    static logging& get()
    {
        return logging_;
    }

private:
    typedef boost::log::sinks::text_ostream_backend console_backend_type;
    typedef boost::log::sinks::synchronous_sink<console_backend_type> console_sink_type;
    typedef boost::log::sinks::text_file_backend file_backend_type;
    typedef boost::log::sinks::synchronous_sink<file_backend_type> file_sink_type;

    /**
     * Opens log to console with level logging::info if compiled
     * without debugging (-DNDEBUG), and logging::debug otherwise.
     *
     * This method is declared as private to ensure that logging
     * remains a singleton instance.
     */
    logging();
    /** set log output format of backend */
    template <typename backend_type>
    void set_formatter(boost::shared_ptr<backend_type> backend);

    /** console log sink */
    boost::shared_ptr<console_sink_type> console_;
    /** file log sink */
    boost::shared_ptr<file_sink_type> file_;
    /** singleton instance */
    static logging logging_;
};

/**
 * Logger source type with severity levels
 *
 * This type may be used by modules to declare there own logger.
 */
typedef boost::log::sources::severity_logger<logging::severity_level> logger;

/**
 * Default logger source
 *
 * We provide a default logger source for use with the LOG* macros.
 * For maximum convenience, this logger is declared as a variable
 * instead of function. If a module should use its own logger, e.g.
 * to add attributes such as the module name, one may conveniently
 * declare a class member logger_, which overrides halmd::logger_.
 *
 * As module loggers will be created in Lua to automatically add
 * attributes such as the module name, we declare logger_ as a
 * shared_ptr, which allows passing it to the module constructor
 * (loggers have no copy constructor). As a safety measure we
 * declare the shared_ptr as const, while the logger itself is
 * mutable.
 */
extern boost::shared_ptr<logger> const logger_;

#define HALMD_LOG(level, format)                        \
{                                                       \
    using namespace halmd;                              \
    BOOST_LOG_SEV(*logger_, level) << format;           \
}                                                       \

#define HALMD_LOG_ONCE(level, format)                   \
{                                                       \
    static bool logged = false;                         \
    if (!logged) {                                      \
        HALMD_LOG(level, format)                        \
        logged = true;                                  \
    }                                                   \
}

#define LOG_FATAL(format)           HALMD_LOG(logging::fatal, format)
#define LOG_FATAL_ONCE(format)      HALMD_LOG_ONCE(logging::fatal, format)
#define LOG_ERROR(format)           HALMD_LOG(logging::error, format)
#define LOG_ERROR_ONCE(format)      HALMD_LOG_ONCE(logging::error, format)
#define LOG_WARNING(format)         HALMD_LOG(logging::warning, format)
#define LOG_WARNING_ONCE(format)    HALMD_LOG_ONCE(logging::warning, format)
#define LOG(format)                 HALMD_LOG(logging::info, format)
#define LOG_ONCE(format)            HALMD_LOG_ONCE(logging::info, format)
#ifndef NDEBUG
# define LOG_DEBUG(format)          HALMD_LOG(logging::debug, format)
# define LOG_DEBUG_ONCE(format)     HALMD_LOG_ONCE(logging::debug, format)
# define LOG_TRACE(format)          HALMD_LOG(logging::trace, format)
# define LOG_TRACE_ONCE(format)     HALMD_LOG_ONCE(logging::trace, format)
#else
# define LOG_DEBUG(format)
# define LOG_DEBUG_ONCE(format)
# define LOG_TRACE(format)
# define LOG_TRACE_ONCE(format)
#endif

} // namespace halmd

#endif /* ! HALMD_IO_LOGGER_HPP */
