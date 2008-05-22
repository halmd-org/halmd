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


// use a cache string to make message formatting faster
BOOST_LOG_FORMAT_MSG(boost::logging::optimize::cache_string_one_str<>)

#define L_(lvl) BOOST_LOG_USE_LOG_IF_LEVEL(::mdsim::log::g_l(), ::mdsim::log::g_l_filter(), lvl)

namespace mdsim { namespace log
{

BOOST_DECLARE_LOG(g_l, boost::logging::logger_format_write<>)
BOOST_DECLARE_LOG_FILTER(g_l_filter, boost::logging::level::holder)

void init(options const& opts);

}} // namespace mdsim::log

#endif /* ! MDSIM_LOG_HPP */
