/* Accumulator with statistical evaluation functions
 *
 * Copyright © 2008-2009  Peter Colberg
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

#ifndef LJGPU_MATH_ACCUM_HPP
#define LJGPU_MATH_ACCUM_HPP

#include <cmath>
#include <limits>
#include <stdint.h>

namespace ljgpu
{

/**
 * Accumulator with statistical evaluation functions
 */
template <typename T>
class accumulator
{
public:
    typedef T value_type;

public:
    accumulator() : n_(0), m_(0), v_(0) {}

    accumulator(uint64_t n, T const& m, T const& v) : n_(n), m_(m), v_(v) {}

    /**
     * accumulate a value
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
	const T t = val - m_;
	n_++;
	m_ += t / n_;
	v_ += t * (val - m_);
	return *this;
    }

    /**
     * accumulate values of another accumulator
     */
    accumulator<T>& operator +=(accumulator<T> const& acc)
    {
	const uint64_t n = n_ + acc.n_;
	v_ = v_ + acc.v_ + std::pow(m_ - acc.m_, 2) * n_ * acc.n_ / n;
	m_ = (n_ * m_ + acc.n_ * acc.m_) / n;
	n_ = n;
	return *this;
    }

    /**
     * reset accumulator to empty state
     */
    void clear()
    {
	n_ = 0;
	m_ = 0;
	v_ = 0;
    }

    /**
     * get accumulator value count
     */
    uint64_t const& count() const
    {
	return n_;
    }

    /**
     * compute mean average
     */
    T mean() const
    {
	if (n_ < 1) {
	    return std::numeric_limits<T>::quiet_NaN();
	}
	return m_;
    }

    /**
     * compute variance
     */
    T var() const
    {
	if (n_ < 2) {
	    return std::numeric_limits<T>::quiet_NaN();
	}
	return v_;
    }

    /**
     * compute standard deviation
     */
    T std() const
    {
	if (n_ < 2) {
	    return std::numeric_limits<T>::quiet_NaN();
	}
	return std::sqrt(v_ / (n_ - 1));
    }

    /**
     * compute standard error of mean
     */
    T err() const
    {
	if (n_ < 2) {
	    return std::numeric_limits<T>::quiet_NaN();
	}
	return sqrt(v_ / (n_ - 1) / n_);
    }

private:
    /** count */
    uint64_t n_;
    /** mean */
    T m_;
    /** variance */
    T v_;
};

} // namespace ljgpu

#endif /* ! LJGPU_MATH_ACCUM_HPP */
