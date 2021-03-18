/*
 * Copyright © 2013       Nicolas Höft
 * Copyright © 2008-2010  Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#ifndef HALMD_NUMERIC_ACCUMULATOR_HPP
#define HALMD_NUMERIC_ACCUMULATOR_HPP

#include <halmd/config.hpp>

#include <boost/type_traits/is_convertible.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/type_traits/has_multiplies.hpp>
#include <boost/utility/enable_if.hpp>
#ifndef __CUDACC__
# include <cmath>
# include <stdint.h>
#endif

namespace halmd {
namespace numeric {
namespace detail {

//
// In order to make the accumulator work with fixed_vector and
// matrices that do not implement operator*(), the multiplication has been
// replaced with element_prod(). For types like double or float that
// do no have an element_prod() implementation, provide one using
// simple multiplication.
//
template <typename T>
inline HALMD_GPU_ENABLED
typename boost::enable_if<boost::has_multiplies<T>, T>::type
element_prod(T const& v, T const& w)
{
    return v * w;
}


/**
 * Accumulator with statistical evaluation functions
 */
template <typename T>
class accumulator
{
public:
    typedef T value_type;
#ifdef __CUDACC__
    typedef size_t size_type; // 64 bit for Fermi in *future* CUDA, 32 bit otherwise
#else
    typedef uint64_t size_type;
#endif

    /**
     * initialize accumulator
     */
    HALMD_GPU_ENABLED accumulator()
      : n_(0), m_(0), v_(0) {}

    /**
     * copy accumulator
     */
    template <typename U>
    HALMD_GPU_ENABLED accumulator(
        accumulator<U> const& acc
      , typename boost::enable_if<boost::is_convertible<U, T> >::type* dummy = 0
    ) : n_(acc.n_), m_(acc.m_), v_(acc.v_) {}

    /**
     * Add value to accumulator.
     */
    template <typename V>
    HALMD_GPU_ENABLED typename boost::enable_if<boost::is_convertible<V, T>, void>::type
    operator()(V const& value)
    {
        //
        // The following method for calculating means and standard
        // deviations with floating point arithmetic is described in
        //
        // D.E. Knuth, Art of Computer Programming, Volume 2: Seminumerical
        // Algorithms, 3rd Edition, 1997, Addison-Wesley, p. 232
        //
        T const t = static_cast<T>(value) - m_;
        n_++;
        m_ += t / n_;
        // see above for use of element_prod()
        v_ += element_prod(t, (static_cast<T>(value) - m_));
    }

    /**
     * Add accumulator to accumulator.
     */
    HALMD_GPU_ENABLED void
    operator()(accumulator<T> const& acc)
    {
        // We have to ensure that the sum of both accumulator counts is
        // non-zero, to avoid division by zero. It suffices to check
        // that the count of the *added* accumulator is non-zero, as
        // it may be discarded if empty.
        if (acc.n_ > 0) {
            typename accumulator<T>::size_type const count = n_ + acc.n_;
            v_ += acc.v_;
            T const diff = m_ - acc.m_;
            v_ += element_prod(diff, diff) * n_ * acc.n_ / count;
            m_ = (n_ * m_ + acc.n_ * acc.m_) / count;
            n_ = count;
        }
    }

    /**
     * Reset accumulator.
     */
    HALMD_GPU_ENABLED void reset()
    {
        *this = accumulator();
    }

    /** count */
    size_type n_;
    /** mean */
    value_type m_;
    /** variance × count */
    value_type v_;
};

/**
 * Returns number of accumulated values.
 */
template <typename T>
HALMD_GPU_ENABLED typename accumulator<T>::size_type const&
count(accumulator<T> const& acc)
{
    return acc.n_;
}

/**
 * Returns mean value.
 */
template <typename T>
HALMD_GPU_ENABLED T const& mean(accumulator<T> const& acc)
{
    return acc.m_;
}

/**
 * Returns variance.
 */
template <typename T>
HALMD_GPU_ENABLED T variance(accumulator<T> const& acc)
{
    return acc.v_ / acc.n_;
}

#ifndef __CUDACC__
using std::sqrt;
#endif

/**
 * compute standard deviation
 */
template <typename T>
HALMD_GPU_ENABLED T sigma(accumulator<T> const& acc)
{
    return sqrt(acc.v_ / acc.n_);
}

/**
 * Returns standard error of mean.
 */
template <typename T>
HALMD_GPU_ENABLED T error_of_mean(accumulator<T> const& acc)
{
    return sqrt((acc.v_ / acc.n_) / (acc.n_ - 1));
}

/**
 * Returns the sum of the accumulated values.
 */
template <typename T>
HALMD_GPU_ENABLED T sum(accumulator<T> const& acc)
{
    return mean(acc) * count(acc);
}

} // namespace detail
} // namespace numeric

// import into top-level namespace
template <typename T>
using accumulator = halmd::numeric::detail::accumulator<T>;

} // namespace halmd

#endif /* ! HALMD_NUMERIC_ACCUMULATOR_HPP */
