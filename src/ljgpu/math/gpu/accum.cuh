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

#ifndef LJGPU_MATH_GPU_ACCUM_CUH
#define LJGPU_MATH_GPU_ACCUM_CUH

namespace ljgpu { namespace gpu
{

/**
 * Accumulator with statistical evaluation functions
 */
template <typename T>
class accumulator
{
public:
    __device__ __host__ inline accumulator() : n_(0), m_(0), v_(0) {}

    __device__ __host__ inline accumulator(unsigned int n, T const& m, T const& v) : n_(n), m_(m), v_(v) {}

    /**
     * accumulate a value
     */
    __device__ __host__ inline accumulator<T>& operator+=(T const& val)
    {
	//
	// The following method for calculating means and standard
	// deviations with floating point arithmetic is described in
	//
	// D.E. Knuth, Art of Computer Programming, Volume 2: Seminumerical
	// Algorithms, 3rd Edition, 1997, Addison-Wesley, p. 232
	//
	T const t = val - m_;
	n_++;
	m_ += t / T(n_);
	v_ += t * (val - m_);
	return *this;
    }

    /**
     * accumulate values of another accumulator
     */
    __device__ __host__ inline accumulator<T>& operator +=(accumulator<T> const& acc)
    {
	unsigned int const n = n_ + acc.n_;
	T const d = m_ - acc.m_;
	v_ = v_ + acc.v_ + d * d * T(n_) * T(acc.n_) / T(n);
	m_ = (T(n_) * m_ + T(acc.n_) * acc.m_) / T(n);
	n_ = n;
	return *this;
    }

    /**
     * get accumulator value count
     */
    __device__ __host__ inline unsigned int count() const
    {
	return n_;
    }

    /**
     * compute mean average
     */
    __device__ __host__ inline T mean() const
    {
	return m_;
    }

    /**
     * compute variance
     */
    __device__ __host__ inline T var() const
    {
	return v_;
    }

private:
    /** count */
    unsigned int n_;
    /** mean */
    T m_;
    /** variance */
    T v_;
};

}} // namespace ljgpu::gpu

#endif /* ! LJGPU_MATH_GPU_ACCUM_CUH */
