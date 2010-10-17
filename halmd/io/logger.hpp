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

#ifndef HALMD_IO_LOGGER_HPP
#define HALMD_IO_LOGGER_HPP

#include <boost/log/sinks/sync_frontend.hpp>
#include <boost/log/sinks/text_file_backend.hpp>
#include <boost/log/sinks/text_ostream_backend.hpp>
#include <boost/log/sources/record_ostream.hpp>
#include <boost/log/sources/severity_logger.hpp>


namespace halmd
{

class logger
{
public:
    enum severity_level
    {
        fatal,
        error,
        warning,
        info,
        debug,
        trace,
    };

    typedef boost::log::sinks::text_ostream_backend console_backend;
    typedef boost::log::sinks::synchronous_sink<console_backend> console_sink;
    typedef boost::log::sinks::text_file_backend file_backend;
    typedef boost::log::sinks::synchronous_sink<file_backend> file_sink;

    logger();
    ~logger();
    void log_to_console(severity_level level);
    void log_to_file(severity_level level, std::string file_name);

    static boost::log::sources::severity_logger<severity_level>& get()
    {
       return logger_;
    }

private:
    static boost::log::sources::severity_logger<severity_level> logger_;

    boost::shared_ptr<console_sink> console_;
    boost::shared_ptr<file_sink> file_;
};

#define __HALMD_LOG__(__level__, __format__)            \
    do {                                                \
        BOOST_LOG_SEV(                                  \
            ::halmd::logger::get()                      \
          , ::halmd::logger::__level__                  \
          ) << __format__;                              \
    } while(0)

#define __HALMD_LOG_ONCE__(__level__, __format__)       \
    do {                                                \
        static bool __logged__ = false;                 \
        if (!__logged__) {                              \
            BOOST_LOG_SEV(                              \
                ::halmd::logger::get()                  \
              , ::halmd::logger::__level__              \
              ) << __format__;                          \
            __logged__ = true;                          \
        }                                               \
    } while(0)

#define LOG_FATAL(__format__)           __HALMD_LOG__(fatal, __format__)
#define LOG_FATAL_ONCE(__format__)      __HALMD_LOG_ONCE__(fatal, __format__)
#define LOG_ERROR(__format__)           __HALMD_LOG__(error, __format__)
#define LOG_ERROR_ONCE(__format__)      __HALMD_LOG_ONCE__(error, __format__)
#define LOG_WARNING(__format__)         __HALMD_LOG__(warning, __format__)
#define LOG_WARNING_ONCE(__format__)    __HALMD_LOG_ONCE__(warning, __format__)
#define LOG(__format__)                 __HALMD_LOG__(info, __format__)
#define LOG_ONCE(__format__)            __HALMD_LOG_ONCE__(info, __format__)

#ifndef NDEBUG
# define LOG_DEBUG(__format__)          __HALMD_LOG__(debug, __format__)
# define LOG_DEBUG_ONCE(__format__)     __HALMD_LOG_ONCE__(debug, __format__)
# define LOG_TRACE(__format__)          __HALMD_LOG__(trace, __format__)
# define LOG_TRACE_ONCE(__format__)     __HALMD_LOG_ONCE__(trace, __format__)
#else
# define LOG_DEBUG(__format__)
# define LOG_DEBUG_ONCE(__format__)
# define LOG_TRACE(__format__)
# define LOG_TRACE_ONCE(__format__)
#endif

} // namespace halmd

#endif /* ! HALMD_IO_LOGGER_HPP */
