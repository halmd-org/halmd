/* Accumulator with statistical evaluation functions
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

#ifndef MDSIM_ACCUMULATOR_HPP
#define MDSIM_ACCUMULATOR_HPP

#include <vector>
#include <math.h>


namespace mdsim
{

/**
 * Accumulator with statistical evaluation functions
 */
template <typename T, typename Allocator = std::allocator<T> >
class accumulator : public std::vector<T, Allocator>
{
public:
    /**
     * accumulate value
     */
    accumulator<T>& operator=(T const& val)
    {
	push_back(val);
	return *this;
    }

    /**
     * compute mean average
     */
    T mean() const
    {
	T val = 0;
	typename std::vector<T, Allocator>::const_iterator it;

	for (it = this->begin(); it != this->end(); it++) {
	    val += *it;
	}
	val /= this->size();

	return val;
    }
};

} // namespace mdsim

#endif /* ! MDSIM_ACCUMULATOR_HPP */
