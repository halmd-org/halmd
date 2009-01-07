/* Statistical evaluation functions
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

#ifndef MDSIM_STATISTICS_HPP
#define MDSIM_STATISTICS_HPP

#include <stddef.h>

namespace mdsim
{

/**
 * compute mean average
 */
template <typename input_iterator>
typename input_iterator::value_type mean(input_iterator const& first, input_iterator const& last)
{
    typename input_iterator::value_type mean = 0;
    size_t count = 0;

    for (input_iterator it = first; it != last; ++it) {
	mean += (*it - mean) / ++count;
    }
    return mean;
}

}

#endif /* ! MDSIM_STATISTICS_HPP */
