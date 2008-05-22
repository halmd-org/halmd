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

BOOST_DEFINE_LOG(g_l, boost::logging::logger_format_write<>)
BOOST_DEFINE_LOG_FILTER(g_l_filter, boost::logging::level::holder)

/**
 * initialize logging
 */
void init(options const& opts) {
    // add formatters and destinations
    g_l()->writer().add_formatter(boost::logging::formatter::high_precision_time("[$dd-$MM-$yyyy $hh:$mm:$ss.$mili] "));
    g_l()->writer().add_formatter(boost::logging::formatter::append_newline());

    // add destinations
    g_l()->writer().add_destination(boost::logging::destination::cout());
    g_l()->writer().add_destination(boost::logging::destination::file(opts.logfile()));
    g_l()->turn_cache_off();

    // set console logging verbosity
    if (opts.verbosity() > 1) {
	g_l_filter()->set_enabled(boost::logging::level::debug);
    }
    else if (opts.verbosity() == 1) {
	g_l_filter()->set_enabled(boost::logging::level::info);
    }
    else {
	g_l_filter()->set_enabled(boost::logging::level::error);
    }
}

}} // namespace mdsim::log
