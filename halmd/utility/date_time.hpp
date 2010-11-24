/*
 * Copyright © 2010  Peter Colberg
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

#ifndef HALMD_UTILITY_DATE_TIME_HPP
#define HALMD_UTILITY_DATE_TIME_HPP

#include <boost/date_time/local_time/local_date_time.hpp>
#include <sstream>
#include <string>

namespace halmd
{

/**
 * returns formatted current local time
 */
inline std::string format_local_time(std::string const& format)
{
    std::ostringstream stream;
    boost::posix_time::ptime const lt = boost::posix_time::second_clock::local_time();
    boost::posix_time::time_facet *const f = new boost::posix_time::time_facet(format.c_str());
    stream.imbue(std::locale(stream.getloc(),f));
    stream << lt;
    return stream.str();
}

} // namespace halmd

#endif /* ! HALMD_UTILITY_DATE_TIME_HPP */
