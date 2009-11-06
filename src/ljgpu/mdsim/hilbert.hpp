/* Lennard-Jones fluid simulation
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

#ifndef LJGPU_MDSIM_HILBERT_HPP
#define LJGPU_MDSIM_HILBERT_HPP

#include <algorithm>
#include <cmath>
#include <ljgpu/math/vector2d.hpp>
#include <ljgpu/math/vector3d.hpp>
#include <ljgpu/util/log.hpp>

namespace ljgpu
{

/**
 * Hilbert space-filling curve generation
 */
template <typename T, unsigned dimension>
class hilbert_sfc
{
private:
    /** periodic box length */
    T const box;
    /** Hilbert space-filling curve recursion depth */
    unsigned int const depth;

public:
    hilbert_sfc(T box, unsigned int depth) : box(box), depth(depth) {}

    /**
     * map 3-/2-dimensional point to 1-dimensional point on Hilbert space curve
     */
    unsigned int operator()(vector<T, dimension> r)
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
        unsigned int const n = 1UL << depth;
        // fractional index of particle's Hilbert cell in [0, n)
        r /= box;
        r = (r - floor(r)) * n;

        // round particle position to center of cell in unit coordinates
        r = (floor(r) + vector<T, dimension>(0.5)) / n;
        // use symmetric coordinates
        r -= vector<T, dimension>(0.5);

        //
        // Jun Wang & Jie Shan, Space-Filling Curve Based Point Clouds Index,
        // GeoComputation, 2005
        //

        // Hilbert code for particle
        unsigned int hcode = 0;

        if (dimension == 3) {
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

#define MASK ((1 << 3) - 1)

            // 32-bit integer for 3D Hilbert code allows a maximum of 10 levels
            for (unsigned int i = 0; i < depth; ++i) {
                // determine Hilbert vertex closest to particle
                vector<unsigned int, dimension> x;
                x[0] = std::signbit(r[0]) & 1;
                x[1] = std::signbit(r[1]) & 1;
                x[2] = std::signbit(r[2]) & 1;
                // lookup Hilbert code
                unsigned int const v = (vc >> (3 * (x[0] + (x[1] << 1) + (x[2] << 2))) & MASK);

                // scale particle coordinates to subcell
                r *= 2;
                r += vector<T, dimension>(x) - 0.5;
                // apply permutation rule according to Hilbert code
                if (v == 0) {
                    vertex_swap(vc, b, h, MASK);
                    vertex_swap(vc, c, e, MASK);
                }
                else if (v == 1 || v == 2) {
                    vertex_swap(vc, c, g, MASK);
                    vertex_swap(vc, d, h, MASK);
                }
                else if (v == 3 || v == 4) {
                    vertex_swap(vc, a, c, MASK);
#ifdef USE_HILBERT_ALT_3D
                    vertex_swap(vc, b, d, MASK);
                    vertex_swap(vc, e, g, MASK);
#endif
                    vertex_swap(vc, f, h, MASK);
                }
                else if (v == 5 || v == 6) {
                    vertex_swap(vc, a, e, MASK);
                    vertex_swap(vc, b, f, MASK);
                }
                else if (v == 7) {
                    vertex_swap(vc, a, g, MASK);
                    vertex_swap(vc, d, f, MASK);
                }

                // add vertex code to partial Hilbert code
                hcode = (hcode << 3) + v;
            }
#undef MASK
        }
        else {
            // Hilbert code-to-vertex lookup table
            unsigned int a = 6;
            unsigned int b = 4;
            unsigned int c = 0;
            unsigned int d = 2;
            // Hilbert vertex-to-code lookup table
            unsigned int vc = 1U << b ^ 2U << c ^ 3U << d;

#define MASK ((1 << 2) - 1)

            // 32-bit integer for 2D Hilbert code allows a maximum of 16 levels
            for (unsigned int i = 0; i < depth; ++i) {
                // determine Hilbert vertex closest to particle
                vector<unsigned int, dimension> x;
                x[0] = std::signbit(r[0]) & 1;
                x[1] = std::signbit(r[1]) & 1;
                // lookup Hilbert code
                unsigned int const v = (vc >> (2 * (x[0] + (x[1] << 1))) & MASK);

                // scale particle coordinates to subcell
                r *= 2;
                r += vector<T, dimension>(x) - 0.5;
                // apply permutation rule according to Hilbert code
                if (v == 0) {
                    vertex_swap(vc, b, d, MASK);
                }
                else if (v == 3) {
                    vertex_swap(vc, a, c, MASK);
                }

                // add vertex code to partial Hilbert code
                hcode = (hcode << 2) + v;
            }
#undef MASK
        }
        return hcode;
    }

private:
    /**
     * swap Hilbert spacing-filling curve vertices
     */
    void vertex_swap(unsigned int& v, unsigned int& a, unsigned int& b, unsigned int mask)
    {
        // swap bits comprising Hilbert codes in vertex-to-code lookup table
        unsigned int const va = ((v >> a) & mask);
        unsigned int const vb = ((v >> b) & mask);
        v = v ^ (va << a) ^ (vb << b) ^ (va << b) ^ (vb << a);
        // update code-to-vertex lookup table
        std::swap(a, b);
    }
};

} // namespace ljgpu

#endif /* ! LJGPU_MDSIM_HILBERT_HPP */
