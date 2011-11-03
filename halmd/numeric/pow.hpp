/*
 * Copyright © 2010  Felix Höfling
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

#ifndef HALMD_NUMERIC_POW_HPP
#define HALMD_NUMERIC_POW_HPP

#include <halmd/config.hpp>

namespace halmd {

/**
 * efficient (but not most) exponentiation with integer power
 * provide exponent as template parameter, adopted from cmath.tcc of libstdc++
 *
 * for algorithm see:
 * D.E. Knuth, Art of Computer Programming, Volume 2: Seminumerical
 * Algorithms, 2nd Edition, 1981, Addison-Wesley, p. 442
 */
template <unsigned int n, typename float_type>
inline HALMD_GPU_ENABLED float_type fixed_pow(float_type x)
{
    float_type y = (n % 2) ? x : 1;

    unsigned _n = n;
    while (_n >>= 1)
    {
        x = x * x;
        if (_n % 2)
            y = y * x;
    }

    return y;
}

/**
  * CUDA-enabled pow() with a priori unknown integer exponent
  */
template <typename float_type>
inline HALMD_GPU_ENABLED float_type pow(float_type x, unsigned short n)
{
    float_type y = (n % 2) ? x : 1;

    while (n >>= 1)
    {
        x = x * x;
        if (n % 2)
            y = y * x;
    }

    return y;
}

} // namespace halmd

#endif /* ! HALMD_NUMERIC_POW_HPP */
