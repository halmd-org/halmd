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

#include <boost/log/sources/record_ostream.hpp>
#include <boost/log/sources/severity_logger.hpp>

#include <halmd/utility/options.hpp>

namespace halmd
{

enum severity_level
{
    trace,
    debug,
    info,
    warning,
    error,
    fatal,
};

extern boost::log::sources::severity_logger<severity_level> logger_;

#define HALMD_LOG(lvl, fmt) \
    do { \
        BOOST_LOG_SEV(logger_, (lvl)) << fmt; \
    } while(0)

#define HALMD_LOG_ONCE(lvl, fmt) \
    do { \
        static bool __logged__ = false; \
        if (!__logged__) { \
            BOOST_LOG_SEV(logger_, (lvl)) << fmt; \
            __logged__ = true; \
        } \
    } while(0)

#define LOG_FATAL(fmt)          HALMD_LOG(fatal, fmt)
#define LOG_FATAL_ONCE(fmt)     HALMD_LOG_ONCE(fatal, fmt)
#define LOG_ERROR(fmt)          HALMD_LOG(error, fmt)
#define LOG_ERROR_ONCE(fmt)     HALMD_LOG_ONCE(error, fmt)
#define LOG_WARNING(fmt)        HALMD_LOG(warning, fmt)
#define LOG_WARNING_ONCE(fmt)   HALMD_LOG_ONCE(warning, fmt)
#define LOG(fmt)                HALMD_LOG(info, fmt)
#define LOG_ONCE(fmt)           HALMD_LOG_ONCE(info, fmt)

#ifndef NDEBUG
# define LOG_DEBUG(fmt)         HALMD_LOG(debug, fmt)
# define LOG_DEBUG_ONCE(fmt)    HALMD_LOG_ONCE(debug, fmt)
# define LOG_TRACE(fmt)         HALMD_LOG(trace, fmt)
# define LOG_TRACE_ONCE(fmt)    HALMD_LOG_ONCE(trace, fmt)
#else
# define LOG_DEBUG(fmt)
# define LOG_DEBUG_ONCE(fmt)
# define LOG_TRACE(fmt)
# define LOG_TRACE_ONCE(fmt)
#endif

namespace io { namespace logger
{

extern void init(po::options const& vm);

}} // namespace io::logger

} // namespace halmd

#endif /* ! HALMD_IO_LOGGER_HPP */
