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

#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>
using namespace boost;

namespace ljgpu { namespace cu
{

/**
 * transformation type
 */
struct accumulate_;

/**
 * accumulate a value
 */
template <typename transform_, typename T>
__device__ typename enable_if<is_same<transform_, accumulate_>, void>::type
transform(unsigned int& n, T& m, T& v, T const& val)
{
    //
    // The following method for calculating means and standard
    // deviations with floating point arithmetic is described in
    //
    // D.E. Knuth, Art of Computer Programming, Volume 2: Seminumerical
    // Algorithms, 3rd Edition, 1997, Addison-Wesley, p. 232
    //
    T const t = val - m;
    n++;
    m += t / static_cast<T>(n);
    v += t * (val - m);
}

/**
 * accumulate values of another accumulator
 */
template <typename transform_, typename T>
__device__ typename enable_if<is_same<transform_, accumulate_>, void>::type
transform(unsigned int& n, T& m, T& v, unsigned int n_, T const& m_, T const& v_)
{
    if (n_ > 0) {
        unsigned int const s = n + n_;
        T const d = m - m_;
        v = v + v_ + d * d * T(n) * T(n_) / T(s);
        m = (T(n) * m + T(n_) * m_) / T(s);
        n = s;
    }
}

}} // namespace ljgpu::cu

#endif /* ! LJGPU_MATH_GPU_ACCUM_CUH */
