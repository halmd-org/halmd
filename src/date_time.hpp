/* Date-time functions
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

#ifndef MDSIM_DATE_TIME_HPP
#define MDSIM_DATE_TIME_HPP

#include <string>
#include <time.h>


namespace mdsim
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

} // namespace mdsim

#endif /* ! MDSIM_DATE_TIME_HPP */
