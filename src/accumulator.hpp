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
    accumulator() : count_(0), m_(0), s_(0)
    {
    }

    /**
     * accumulate single value
     */
    accumulator<T>& operator+=(T const& val)
    {
	//
	// The following method for calculating means and standard
	// deviations with floating point arithmetic is described in
	//
	// D.E. Knuth, Art of Computer Programming, Volume 2: Seminumerical
	// Algorithms, 3rd Edition, 1997, Addison-Wesley, p. 232
	//

	++count_;

	T t = val - m_;
	m_ += t / count_;
	s_ += t * (val - m_);

	return *this;
    }

    /**
     * reset accumulator to empty state
     */
    void clear()
    {
	count_ = 0;
	m_ = 0;
	s_ = 0;
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
	assert(count_ > 0);
	return m_;
    }

    /**
     * compute standard deviation
     */
    T std() const
    {
	assert(count_ > 1);
	return sqrt(s_ / (count_ - 1.));
    }

    /**
     * compute standard error of mean
     */
    T err() const
    {
	assert(count_ > 1);
	return sqrt(s_ / (count_ - 1.) / count_);
    }

private:
    size_t count_;
    T m_, s_;
};

} // namespace mdsim

#endif /* ! MDSIM_ACCUMULATOR_HPP */
