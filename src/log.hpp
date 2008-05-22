/* Logging
 *
 * Copyright (C) 2008  Peter Colberg
 *
 * This program is free software: you can redistribute it and/or modify
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

#ifndef MDSIM_LOG_HPP
#define MDSIM_LOG_HPP

#include <boost/logging/format_fwd.hpp>
#include "options.hpp"


/** log informational messages */
#define LOG(fmt) L_(logger) << fmt

/** log error messages */
#define LOG_ERROR(fmt) L_(logger_error) << "[ERROR] " << fmt

/** log debug-level messages */
#ifndef NDEBUG
# define LOG_DEBUG(fmt) L_(logger_debug) << "[DEBUG] " << fmt
# define SCOPED_LOG_DEBUG(fmt) BOOST_SCOPED_LOG_CTX(L_(logger_debug)) << "[DEBUG] " << fmt
#else
# define LOG_DEBUG(fmt)
# define SCOPED_LOG_DEBUG(fmt)
#endif


// use a cache string to make message formatting faster
BOOST_LOG_FORMAT_MSG(boost::logging::optimize::cache_string_one_str<>)

namespace mdsim { namespace log
{

void init(options const& opts);

typedef boost::logging::scenario::ts::use<
    boost::logging::scenario::ts::filter_::none,
    boost::logging::scenario::ts::level_::no_levels,
    boost::logging::scenario::ts::logger_::none> finder;

BOOST_DECLARE_LOG_FILTER(log_filter, finder::filter)

BOOST_DECLARE_LOG(logger, finder::logger)
BOOST_DECLARE_LOG(logger_error, finder::logger)
#ifndef NDEBUG
BOOST_DECLARE_LOG(logger_debug, finder::logger)
#endif

}} // namespace mdsim::log

#define L_(logger) BOOST_LOG_USE_LOG_IF_FILTER(::mdsim::log::logger(), ::mdsim::log::log_filter()->is_enabled())

#endif /* ! MDSIM_LOG_HPP */
