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

#include <boost/logging/format.hpp>
#include <boost/logging/format/formatter/high_precision_time.hpp>
#include "log.hpp"


namespace mdsim { namespace log
{

/**
 * initialize logging
 */
void init(options const& opts) {
    // use microsecond-resolution log timestamps
    boost::logging::formatter::high_precision_time hpt("[$dd-$MM-$yyyy $hh:$mm:$ss.$micro] ");

    // add log formatters
    logger()->writer().add_formatter(hpt);
    logger()->writer().add_formatter(boost::logging::formatter::append_newline());
    logger_error()->writer().add_formatter(hpt);
    logger_error()->writer().add_formatter(boost::logging::formatter::append_newline());
#ifndef NDEBUG
    logger_debug()->writer().add_formatter(hpt);
    logger_debug()->writer().add_formatter(boost::logging::formatter::append_newline());
#endif

    boost::logging::destination::file logfile(opts.logfile());

    // output informational messages to file
    logger()->writer().add_destination(logfile);
    if (opts.verbosity() > 0) {
	// output informational messages to console
	logger()->writer().add_destination(boost::logging::destination::cerr());
    }
    logger()->mark_as_initialized();

    // output error messages to console and file
    logger_error()->writer().add_destination(boost::logging::destination::cerr());
    logger_error()->writer().add_destination(logfile);
    logger_error()->mark_as_initialized();

#ifndef NDEBUG
    if (opts.verbosity() > 1) {
	// output debug-level messages to console and file
	logger_debug()->writer().add_destination(boost::logging::destination::cerr());
	logger_debug()->writer().add_destination(logfile);
    }
    logger_debug()->mark_as_initialized();
#endif
}

BOOST_DEFINE_LOG_FILTER(log_filter, finder::filter)

BOOST_DEFINE_LOG(logger, finder::logger)
BOOST_DEFINE_LOG(logger_error, finder::logger)
#ifndef NDEBUG
BOOST_DEFINE_LOG(logger_debug, finder::logger)
#endif

}} // namespace mdsim::log
