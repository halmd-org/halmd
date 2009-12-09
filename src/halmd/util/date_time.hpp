/* Date-time functions
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

#ifndef HALMD_UTIL_DATE_TIME_HPP
#define HALMD_UTIL_DATE_TIME_HPP

#include <string>
#include <time.h>

namespace halmd
{

/**
 * Date-time functions
 */
class date_time
{
public:
    /**
     * returns formatted current local time
     */
    static std::string format(std::string const& fmt)
    {
        char str[256];
        time_t t;

        // time in seconds since the epoch
        time(&t);
        // format local time
        strftime(str, sizeof(str), fmt.c_str(), localtime(&t));

        return str;
    }
};

} // namespace halmd

#endif /* ! HALMD_UTIL_DATE_TIME_HPP */
