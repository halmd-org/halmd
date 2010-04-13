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

#include <algorithm>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <cmath>

#include <halmd/mdsim/host/sort/hilbert.hpp>
#include <halmd/util/logger.hpp>

using namespace boost;
using namespace std;

namespace halmd { namespace mdsim { namespace host { namespace sort
{

template <int dimension, typename float_type>
hilbert<dimension, float_type>::hilbert(options const& vm)
    : _Base(vm)
    // dependency injection
    , particle(dynamic_pointer_cast<particle_type>(module<mdsim::particle<dimension> >::fetch(vm)))
    , box(dynamic_pointer_cast<box_type>(module<mdsim::box<dimension> >::fetch(vm)))
    , neighbor(dynamic_pointer_cast<neighbor_type>(module<mdsim::neighbor<dimension> >::fetch(vm)))
{
    // set Hilbert space-filling curve recursion depth
    unsigned int ncell = *max_element(neighbor->ncell_.begin(), neighbor->ncell_.end());
    unsigned int depth = static_cast<int>(std::ceil(std::log(static_cast<double>(ncell)) / M_LN2));
    // 32-bit integer for 2D Hilbert code allows a maximum of 16/10 levels
    depth = min((dimension == 3) ? 10U : 16U, depth);

    LOG("Hilbert vertex recursion depth: " << depth);

    // generate 1-dimensional Hilbert curve mapping of cell lists
    typedef std::pair<cell_list*, unsigned int> pair;
    std::vector<pair> pairs;
    ::vector<unsigned int, dimension> x;
    for (x[0] = 0; x[0] < neighbor->ncell_[0]; ++x[0]) {
        for (x[1] = 0; x[1] < neighbor->ncell_[1]; ++x[1]) {
            if (dimension == 3) {
                for (x[2] = 0; x[2] < neighbor->ncell_[2]; ++x[2]) {
                    ::vector<float_type, dimension> r(x);
                    for (size_t i = 0; i < dimension; ++i) {
                        r[i] = (r[i] + 0.5) * neighbor->cell_length_[i];
                    }
                    pairs.push_back(std::make_pair(&neighbor->cell_(x), map(r, depth)));
                }
            }
            else {
                ::vector<float_type, dimension> r(x);
                for (size_t i = 0; i < dimension; ++i) {
                    r[i] = (r[i] + 0.5) * neighbor->cell_length_[i];
                }
                pairs.push_back(std::make_pair(&neighbor->cell_(x), map(r, depth)));
            }
        }
    }
    stable_sort(pairs.begin(), pairs.end(), bind(&pair::second, _1) < bind(&pair::second, _2));
    cell_.clear();
    cell_.reserve(neighbor->cell_.size());
    transform(pairs.begin(), pairs.end(), back_inserter(cell_), bind(&pair::first, _1));
}

/**
 * Order particles after Hilbert space-filling curve
 */
template <int dimension, typename float_type>
void hilbert<dimension, float_type>::order()
{
    // particle binning
    neighbor->update_cells();
    // generate index sequence according to Hilbert-sorted cells
    std::vector<unsigned int> index;
    index.reserve(particle->nbox);
    BOOST_FOREACH(cell_list* cell, cell_) {
        BOOST_FOREACH(unsigned int p, *cell) {
            index.push_back(p);
        }
    }
    // reorder particles in memory
    particle->rearrange(index);
}

/**
 * Map 3-/2-dimensional point to 1-dimensional point on Hilbert space curve
 */
template <int dimension, typename float_type>
unsigned int hilbert<dimension, float_type>::map(vector_type r, unsigned int depth)
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
    float_type const n = 1UL << depth;
    // fractional index of particle's Hilbert cell in [0, n)
    for (size_t i = 0; i < dimension; ++i) {
        r[i] /= box->length()[i];
    }
    r = (r - floor(r)) * n;

    // round particle position to center of cell in unit coordinates
    r = (floor(r) + vector_type(0.5)) / n;
    // use symmetric coordinates
    r -= vector_type(0.5);

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
            ::vector<unsigned int, dimension> x;
            x[0] = std::signbit(r[0]) & 1;
            x[1] = std::signbit(r[1]) & 1;
            x[2] = std::signbit(r[2]) & 1;
            // lookup Hilbert code
            unsigned int const v = (vc >> (3 * (x[0] + (x[1] << 1) + (x[2] << 2))) & MASK);

            // scale particle coordinates to subcell
            r *= 2;
            r += vector_type(x) - 0.5;
            // apply permutation rule according to Hilbert code
            if (v == 0) {
                swap(vc, b, h, MASK);
                swap(vc, c, e, MASK);
            }
            else if (v == 1 || v == 2) {
                swap(vc, c, g, MASK);
                swap(vc, d, h, MASK);
            }
            else if (v == 3 || v == 4) {
                swap(vc, a, c, MASK);
#ifdef USE_HILBERT_ALT_3D
                swap(vc, b, d, MASK);
                swap(vc, e, g, MASK);
#endif
                swap(vc, f, h, MASK);
            }
            else if (v == 5 || v == 6) {
                swap(vc, a, e, MASK);
                swap(vc, b, f, MASK);
            }
            else if (v == 7) {
                swap(vc, a, g, MASK);
                swap(vc, d, f, MASK);
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
            ::vector<unsigned int, dimension> x;
            x[0] = std::signbit(r[0]) & 1;
            x[1] = std::signbit(r[1]) & 1;
            // lookup Hilbert code
            unsigned int const v = (vc >> (2 * (x[0] + (x[1] << 1))) & MASK);

            // scale particle coordinates to subcell
            r *= 2;
            r += vector_type(x) - 0.5;
            // apply permutation rule according to Hilbert code
            if (v == 0) {
                swap(vc, b, d, MASK);
            }
            else if (v == 3) {
                swap(vc, a, c, MASK);
            }

            // add vertex code to partial Hilbert code
            hcode = (hcode << 2) + v;
        }
#undef MASK
    }
    return hcode;
}

/**
 * Swap Hilbert spacing-filling curve vertices
 */
template <int dimension, typename float_type>
void hilbert<dimension, float_type>::swap(unsigned int& v, unsigned int& a, unsigned int& b, unsigned int mask)
{
    // swap bits comprising Hilbert codes in vertex-to-code lookup table
    unsigned int const va = ((v >> a) & mask);
    unsigned int const vb = ((v >> b) & mask);
    v = v ^ (va << a) ^ (vb << b) ^ (va << b) ^ (vb << a);
    // update code-to-vertex lookup table
    std::swap(a, b);
}

// explicit instantiation
#ifndef USE_HOST_SINGLE_PRECISION
template class hilbert<3, double>;
template class hilbert<2, double>;
#else
template class hilbert<3, float>;
template class hilbert<2, float>;
#endif

}}}} // namespace halmd::mdsim::host::sort
