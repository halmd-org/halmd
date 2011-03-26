/*
 * Copyright © 2008-2010  Peter Colberg
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

#ifndef HALMD_MDSIM_SORTS_HILBERT_KERNEL_HPP
#define HALMD_MDSIM_SORTS_HILBERT_KERNEL_HPP

#include <boost/math/special_functions/sign.hpp> // std::signbit is not portable
#include <boost/mpl/equal_to.hpp>
#include <boost/mpl/int.hpp>
#include <boost/utility/enable_if.hpp>
#ifdef __CUDACC__
# include <cuda_runtime.h> // signbit
#else
# include <cmath>
# include <algorithm>
#endif

#include <halmd/config.hpp>
#ifdef __CUDACC__
# include <halmd/algorithm/gpu/bits.cuh> // swap
#endif

namespace halmd
{
namespace mdsim { namespace sorts
{
namespace hilbert_kernel
{
namespace detail
{

// import into detail namespace
#ifdef __CUDACC__
using algorithm::gpu::swap;
#else
using boost::math::signbit;
using std::swap;
#endif

/**
 * Swap Hilbert spacing-filling curve vertices
 */
inline HALMD_GPU_ENABLED void swap(
    unsigned int& v
  , unsigned int& a
  , unsigned int& b
  , unsigned int mask)
{
    // swap bits comprising Hilbert codes in vertex-to-code lookup table
    unsigned int const va = ((v >> a) & mask);
    unsigned int const vb = ((v >> b) & mask);
    v = v ^ (va << a) ^ (vb << b) ^ (va << b) ^ (vb << a);
    // update code-to-vertex lookup table
    swap(a, b);
}

/**
 * Map 3-dimensional point to 1-dimensional point on Hilbert space curve
 */
template <typename vector_type>
HALMD_GPU_ENABLED
typename boost::enable_if<
    boost::mpl::equal_to<
        boost::mpl::int_<vector_type::static_size>
      , boost::mpl::int_<3>
    >
  , unsigned int>::type map(vector_type r, unsigned int depth)
{
    //
    // Jun Wang & Jie Shan, Space-Filling Curve Based Point Clouds Index,
    // GeoComputation, 2005
    //

    // Hilbert code for particle
    unsigned int hcode = 0;

    // Hilbert code-to-vertex lookup table
    unsigned int a = 21;
    unsigned int b = 18;
    unsigned int c = 12;
    unsigned int d = 15;
    unsigned int e = 3;
    unsigned int f = 0;
    unsigned int g = 6;
    unsigned int h = 9;
    // Hilbert vertex-to-code lookup table
    unsigned int vc = 1U << b ^ 2U << c ^ 3U << d ^ 4U << e ^ 5U << f ^ 6U << g ^ 7U << h;

    unsigned int const mask = (1 << 3) - 1;

    // 32-bit integer for 3D Hilbert code allows a maximum of 10 levels
    for (unsigned int i = 0; i < depth; ++i) {
        // determine Hilbert vertex closest to particle
        unsigned int const x = signbit(r[0]) & 1;
        unsigned int const y = signbit(r[1]) & 1;
        unsigned int const z = signbit(r[2]) & 1;
        // lookup Hilbert code
        unsigned int const v = (vc >> (3 * (x + (y << 1) + (z << 2))) & mask);

        // scale particle coordinates to subcell
        r *= 2;
        r[0] += x;
        r[0] -= 0.5f; // promoted to double by compiler, if needed
        r[1] += y;
        r[1] -= 0.5f;
        r[2] += z;
        r[2] -= 0.5f;
        // apply permutation rule according to Hilbert code
        if (v == 0) {
            swap(vc, b, h, mask);
            swap(vc, c, e, mask);
        }
        else if (v == 1 || v == 2) {
            swap(vc, c, g, mask);
            swap(vc, d, h, mask);
        }
        else if (v == 3 || v == 4) {
            swap(vc, a, c, mask);
#ifdef USE_HILBERT_ALT_3D
            swap(vc, b, d, mask);
            swap(vc, e, g, mask);
#endif
            swap(vc, f, h, mask);
        }
        else if (v == 5 || v == 6) {
            swap(vc, a, e, mask);
            swap(vc, b, f, mask);
        }
        else if (v == 7) {
            swap(vc, a, g, mask);
            swap(vc, d, f, mask);
        }

        // add vertex code to partial Hilbert code
        hcode = (hcode << 3) + v;
    }

    return hcode;
}

/**
 * Map 2-dimensional point to 1-dimensional point on Hilbert space curve
 */
template <typename vector_type>
HALMD_GPU_ENABLED
typename boost::enable_if<
    boost::mpl::equal_to<
        boost::mpl::int_<vector_type::static_size>
      , boost::mpl::int_<2>
    >
  , unsigned int>::type map(vector_type r, unsigned int depth)
{
    //
    // Jun Wang & Jie Shan, Space-Filling Curve Based Point Clouds Index,
    // GeoComputation, 2005
    //

    // Hilbert code for particle
    unsigned int hcode = 0;

    // Hilbert code-to-vertex lookup table
    unsigned int a = 6;
    unsigned int b = 4;
    unsigned int c = 0;
    unsigned int d = 2;
    // Hilbert vertex-to-code lookup table
    unsigned int vc = 1U << b ^ 2U << c ^ 3U << d;

    unsigned int const mask = (1 << 2) - 1;

    // 32-bit integer for 2D Hilbert code allows a maximum of 16 levels
    for (unsigned int i = 0; i < depth; ++i) {
        // determine Hilbert vertex closest to particle
        unsigned int const x = signbit(r[0]) & 1;
        unsigned int const y = signbit(r[1]) & 1;
        // lookup Hilbert code
        unsigned int const v = (vc >> (2 * (x + (y << 1))) & mask);

        // scale particle coordinates to subcell
        r *= 2;
        r[0] += x;
        r[0] -= 0.5f; // promoted to double by compiler, if needed
        r[1] += y;
        r[1] -= 0.5f;
        // apply permutation rule according to Hilbert code
        if (v == 0) {
            swap(vc, b, d, mask);
        }
        else if (v == 3) {
            swap(vc, a, c, mask);
        }

        // add vertex code to partial Hilbert code
        hcode = (hcode << 2) + v;
    }

    return hcode;
}

} // namespace detail

/**
 * Map 3-/2-dimensional point to 1-dimensional point on Hilbert space curve
 */
template <typename vector_type>
HALMD_GPU_ENABLED unsigned int map(vector_type r, unsigned int depth)
{
    //
    // We need to avoid ambiguities during the assignment of a particle
    // to a subcell, i.e. the particle position should never lie on an
    // edge or corner of multiple subcells, or the algorithm will have
    // trouble converging to a definite Hilbert curve.
    //
    // Therefore, we use a simple cubic lattice of predefined dimensions
    // according to the number of cells at the deepest recursion level,
    // and round the particle position to the nearest center of a cell.
    //

    // Hilbert cells per dimension at deepest recursion level
    typename vector_type::value_type const n = 1UL << depth;
    // fractional index of particle's Hilbert cell in [0, n)
    r = (r - floor(r)) * n;

    // round particle position to center of cell in unit coordinates
    r = (floor(r) + vector_type(0.5f)) / n;
    // use symmetric coordinates
    r -= vector_type(0.5f);

    return detail::map(r, depth);
}

} // namespace hilbert_kernel

}} // namespace mdsim::sorts

} // namespace halmd

#endif /* ! HALMD_MDSIM_SORTS_HILBERT_KERNEL_HPP */
