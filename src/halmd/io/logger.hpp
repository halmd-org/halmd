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

#include <halmd/options.hpp>

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

#define _LOG(lvl)          BOOST_LOG_SEV(logger_, (lvl))

#ifndef NDEBUG
# define LOG_TRACE(fmt)    _LOG(trace) << fmt
# define LOG_DEBUG(fmt)    _LOG(debug) << fmt
#else
# define LOG_TRACE(fmt)
# define LOG_DEBUG(fmt)
#endif
#define LOG(fmt)           _LOG(info) << fmt
#define LOG_WARNING(fmt)   _LOG(warning) << fmt
#define LOG_ERROR(fmt)     _LOG(error) << fmt
#define LOG_FATAL(fmt)     _LOG(fatal) << fmt

namespace io { namespace logger
{

extern void init(options const& vm);

}} // namespace io::logger

} // namespace halmd

#endif /* ! HALMD_IO_LOGGER_HPP */
