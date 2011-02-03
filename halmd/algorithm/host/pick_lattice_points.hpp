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

#ifndef HALMD_ALGORITHM_HOST_PICK_LATTICE_POINTS_HPP
#define HALMD_ALGORITHM_HOST_PICK_LATTICE_POINTS_HPP

#include <algorithm>
#include <boost/foreach.hpp>
#include <cmath>
#include <map>
#include <vector>

#include <halmd/io/logger.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/numeric/gcd.hpp>

namespace halmd { namespace algorithm { namespace host
{

/**
 *  pick (almost) equally distributed lattice points from a concentric shell
 *
 * @param radii:        list of shell radii
 * @param unit_cell:    edge lengths @f$ (a_1, a_2, ...) @f$ of cuboid unit cell
 * @param tolerance:    relative tolerance on radii, defines thickness of shells
 * @param max_count:    maximal number of points in each shell
 * @returns list of lattice point shells, the points of each shell have the same
 *          distance to the origin (within the given tolerance)
 *
 * The constructed lattice points are of the form @f$ \vec r = n (h a_1, k a_2, l a_3, ...) @f$
 * where the Miller indices @f$ (h, k, l) @f$ are coprime and @f$ n @f$ is integer. Lattice
 * points with a smaller sum @f$ h+k+l @f$ are preferred.
 *
 * Caveat: since the number of points per shell is sharply limited, the returned points
 * are not completely symmetric in space (some possible combinations @f$ (h,k,l) @f$ may
 * be missing for the largest sum @f$ h+k+l @f$)
 */
template <typename float_type, typename vector_type>
std::multimap<float_type, vector_type> pick_lattice_points_from_shell(
    std::vector<float_type> const& radii
  , vector_type const& unit_cell
  , float_type tolerance
  , unsigned int max_count
) // TODO use the InputIterator/OutputIterator idiom
{
    using namespace std;

    enum { dimension = vector_type::static_size };

    typedef fixed_vector<unsigned int, dimension> index_type;
    typedef multimap<float_type, vector_type> map_type;
    map_type lattice_points;

    // determine maximum Miller index that fits in the range of radii
    // (e_i) × (h,k,l) ≤ r_max ⇒ (h,k,l) ≤ r_max / max(e_i)
    float_type r_max = *max_element(radii.begin(), radii.end()) * (1 + tolerance);
    unsigned int miller_max = static_cast<unsigned int>(floor(
        r_max / *max_element(unit_cell.begin(), unit_cell.end())
    ));
    LOG_DEBUG("generate lattice points with a maximum sum of Miller indices: " << miller_max);

    // generate all lattice points bounded h+k+l ≤ miller_max,
    // loop over tuple index (hkl), or (hk) in two dimensions,
    // using idx[2] = h+k+l, idx[1] = h+k, and idx[0] = h as loop variables
    index_type idx(0);
    while (idx[dimension-1] <= miller_max) {
        // construct (hkl) from partial index sums
        index_type hkl;
        hkl[0] = idx[0];
        for (unsigned int j = 1; j < dimension; ++j) {
            hkl[j] = idx[j] - idx[j-1];
        }
        // condition: tuple elements are coprime
        if (is_coprime(hkl))
        {
            // construct and store lattice points along the direction (hkl)
            // with magnitude close to the desired values
            vector_type r0 = element_prod(unit_cell, static_cast<vector_type>(hkl));
            float_type r0_norm = norm_2(r0);
            BOOST_FOREACH (float_type r, radii) {
                if (lattice_points.count(r) < max_count) {
                    // find integer n such that abs(norm_2(n * r0) - r) / r < tolerance
                    // 1) round to nearest integer
                    unsigned int n = floor(r / r0_norm + float_type(.5));
                    // 2) check if this is good enough
                    if (n > 0 && abs(n * r0_norm - r) < r * tolerance) {
                        lattice_points.insert(make_pair(r, n * r0));
                    }
                }
            }
        }
        // increment index tuple at end of loop,
        // obey 0 ≤ idx[0] ≤ idx[1] ≤ ... ≤ miller_max
        ++idx[0];                            // always increment first 'digit' (element)
        for (unsigned int j = 0; j < dimension - 1; ++j) {
            if (idx[j] <= idx[j+1]) {        // increment 'digit' and test upper bound
                break;                       // still within range, exit loop
            }
            idx[j] = 0;                      // reset this 'digit'
            ++idx[j+1];                      // increment next 'digit'
        }
    };
    LOG_TRACE("lattice points constructed: " << lattice_points.size());

    // output list of lattice points for debug trace
    BOOST_FOREACH (float_type r, radii) {
        unsigned int count = lattice_points.count(r);
        LOG_TRACE(count << " lattice points with r ≈ " << r);
        if (!count) {
            LOG_WARNING("No lattice points compatible with r ≈ " << r);
        }
#ifndef NDEBUG
        typedef pair<typename map_type::const_iterator, typename map_type::const_iterator> range_type;
        for (range_type shell = lattice_points.equal_range(r); shell.first != shell.second; ++shell.first) {
            vector_type const& point = shell.first->second;
            index_type hkl = static_cast<index_type>(round(element_div(point, unit_cell)));
            hkl /= greatest_common_divisor(hkl);
            LOG_TRACE("  r = " << norm_2(point) << ", (hkl) = " << hkl);
        }
#endif
    }

    return lattice_points;
}

}}} // namespace halmd::algorithm::host

#endif /* ! HALMD_ALGORITHM_HOST_PICK_LATTICE_POINTS_HPP */
