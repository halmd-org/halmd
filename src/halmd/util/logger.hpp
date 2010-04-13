/* Logging
 *
 * Copyright © 2008-2009  Peter Colberg
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

#ifndef HALMD_UTIL_LOGGER_HPP
#define HALMD_UTIL_LOGGER_HPP

#include <boost/log/sources/record_ostream.hpp>
#include <boost/log/sources/severity_logger.hpp>

namespace halmd { namespace logger
{

enum severity_level
{
    debug,
    info,
    warning,
    error,
};

extern boost::log::sources::severity_logger<severity_level> slg;
extern void init(std::string const& filename, int verbosity);

}} // namespace halmd::logger

#define _LOGGER(level) BOOST_LOG_SEV(halmd::logger::slg, halmd::logger::level)

/** log informational messages */
#define LOG(fmt) _LOGGER(info) << fmt

/** log warning messages */
#define LOG_WARNING(fmt) _LOGGER(warning) << "[WARNING] " << fmt

/** log error messages */
#define LOG_ERROR(fmt) _LOGGER(error) << "[ERROR] " << fmt

/** log debug-level messages */
#ifndef NDEBUG
# define LOG_DEBUG(fmt) _LOGGER(debug) << "[DEBUG] " << fmt
#else
# define LOG_DEBUG(fmt)
#endif

#endif /* ! HALMD_UTIL_LOGGER_HPP */
