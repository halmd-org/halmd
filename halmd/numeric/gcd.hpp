/*
 * Copyright © 2011  Felix Höfling
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

#ifndef HALMD_NUMERIC_GCD_HPP
#define HALMD_NUMERIC_GCD_HPP

#include <boost/array.hpp>

#include <halmd/config.hpp>

namespace halmd
{

/**
 * determine greatest common divisor of two integer numbers using
 * Euclid's algorithm
 *
 * for algorithm see:
 * Wikipedia http://en.wikipedia.org/wiki/Euclidean_algorithm
 *
 * or
 *
 * D.E. Knuth, Art of Computer Programming, Volume 1: Fundamental
 * Algorithms, 3nd ed., 1997, Addison-Wesley.
 */
template <typename T>
inline HALMD_GPU_ENABLED T greatest_common_divisor(T a, T b)
{
    while (b != 0) {
        T t = b;
        b = a % b;
        a = t;
    }
    return a;
}

/** returns true if 'a' and 'b' are coprime numbers */
template <typename T>
inline HALMD_GPU_ENABLED bool is_coprime(T const& a, T const& b)
{
    return greatest_common_divisor(a, b) == 1;
}

//
// overloaded functions for array of numbers
//
/** returns greatest common divisor of all numbers in array */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED T greatest_common_divisor(boost::array<T, N> const& a);

template <typename T>
inline HALMD_GPU_ENABLED T greatest_common_divisor(boost::array<T, 1> const& a)
{
    return a[0];
}

template <typename T, size_t N>
inline HALMD_GPU_ENABLED T greatest_common_divisor(boost::array<T, N> const& a)
{
    // compute gcd(a[0], gcd(a[1], gcd(a[2], ...)))
    boost::array<T, N-1> b;
    std::copy(a.begin() + 1, a.end(), b.begin());

    return greatest_common_divisor(a[0], greatest_common_divisor(b));
}

/** returns true if all numbers in array are coprime */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED bool is_coprime(boost::array<T, N> const& a)
{
    return greatest_common_divisor(a) == 1;
}

} // namespace halmd

#endif /* ! HALMD_NUMERIC_GCD_HPP */
