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

#include <math.h>
#include <assert.h>


namespace mdsim
{

/**
 * Accumulator with statistical evaluation functions
 */
template <typename T>
class accumulator
{
public:
    accumulator() : count_(0), mean_(0), std_(0)
    {
    }

    /**
     * accumulate value
     */
    accumulator<T>& add(T const& val)
    {
	// accumulate mean average
	mean_ += val;
	// accumulate standard deviation
	std_ += val * val;

	count_++;
	return *this;
    }

    /**
     * reset accumulator to empty state
     */
    void clear()
    {
	mean_ = 0;
	std_ = 0;
	count_ = 0;
    }

    /**
     * get accumulator value count
     */
    size_t count() const
    {
	return count_;
    }

    /**
     * compute mean average
     */
    T mean() const
    {
	assert(count_ != 0);
	return mean_ / count_;
    }

    /**
     * compute standard deviation
     */
    T std() const
    {
	assert(count_ != 0);
	return sqrt((std_ - mean_ * mean_ / count_) / count_);
    }

private:
    size_t count_;
    T mean_, std_;
};

} // namespace mdsim

#endif /* ! MDSIM_ACCUMULATOR_HPP */
