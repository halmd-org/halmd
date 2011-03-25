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
#include <boost/concept_check.hpp>
#include <boost/concept/requires.hpp>
#include <boost/type_traits/is_same.hpp>
#include <functional>
#include <numeric>
#include <vector>
#include <utility>

#include <halmd/io/logger.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/numeric/gcd.hpp>

namespace halmd
{
namespace algorithm { namespace host
{

/**
 *  pick (almost) equally distributed lattice points from a concentric shell
 *
 * @param radius_begin: iterator at the begin of a range of shell radii
 * @param radius_end:   iterator at the end of the range
 * @param result:       output iterator over pairs (radius, lattice point)
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
template <typename float_type, typename vector_type, typename InputIterator, typename OutputIterator>
BOOST_CONCEPT_REQUIRES(
    ((boost::InputIterator<InputIterator>))
    ((boost::OutputIterator<OutputIterator, std::pair<float_type const, vector_type> >))
  , (void)) // return type
pick_lattice_points_from_shell(
    InputIterator radius_begin, InputIterator radius_end
  , OutputIterator result
  , vector_type const& unit_cell
  , float_type tolerance
  , unsigned int max_count
)
{
    using namespace std;

    enum { dimension = vector_type::static_size };
    typedef fixed_vector<unsigned int, dimension> index_type;

    // return on empty range
    if (radius_begin == radius_end) return;

    // keep track of the number of constructed lattice points
    vector<unsigned int> count(radius_end - radius_begin, 0u);

    // determine maximum sum of Miller indices that fits in the range of radii
    // (e_i) × (h,k,l) ≤ r_max ⇒ (h,k,l)_i ≤ r_max / e_i ⇒ h+k+l ≤ r_max × sum(1/e_i)
    float_type r_max = *max_element(radius_begin, radius_end) * (1 + tolerance);
    index_type miller_sum_max = static_cast<index_type>(element_div(vector_type(r_max), unit_cell));
    // compute partial sum: miller_sum_max = (h_max, h_max + k_max, h_max + k_max + l_max)
    partial_sum(miller_sum_max.begin(), miller_sum_max.end(), miller_sum_max.begin());
    LOG_DEBUG("generate lattice points with a maximum sum of Miller indices: " << miller_sum_max[dimension-1]);

    // generate all lattice points bounded h+k+l ≤ miller_sum_max,
    // loop over tuple index (hkl), or (hk) in two dimensions,
    // using idx[2] = h+k+l, idx[1] = h+k, and idx[0] = h as loop variables
    index_type idx(0);
    while (idx[dimension-1] <= miller_sum_max[dimension-1]) {
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
            unsigned int i = 0;
            for (InputIterator r_it = radius_begin; r_it != radius_end; ++r_it, ++i) {
                if (count[i] < max_count) {
                    // find integer n such that abs(norm_2(n * r0) - r) / r ≤ tolerance
                    // 1) round to nearest integer
                    unsigned int n = floor(*r_it / r0_norm + float_type(.5));
                    // 2) check if this is good enough
                    if (n > 0 && (abs(n * r0_norm - *r_it) <= *r_it * tolerance)) {
                        vector_type point = n * r0;
                        *result++ = make_pair(*r_it, point);
                        ++count[i];
#ifndef NDEBUG
                        index_type hkl_reduced = hkl / greatest_common_divisor(hkl);
                        LOG_TRACE(
                            "  r = " << norm_2(point) << ", (hkl) = " << hkl_reduced
                         << " ⇒ " << n << " × √" << inner_prod(hkl, hkl)
                        );
#endif
                    }
                }
            }
        }
        // increment index tuple at end of loop,
        // obey 0 ≤ idx[j] ≤ miller_sum_max[j]
        ++idx[0];                            // always increment first 'digit' (element)
        for (unsigned int j = 0; j < dimension - 1; ++j) {
            if (idx[j] <= miller_sum_max[j]) {  // test upper bound
                break;                          // still within range, exit loop
            }
            idx[j] = 0;                         // reset this 'digit'
            ++idx[j+1];                         // increment next 'digit'
        }
    };

    LOG_TRACE("lattice points constructed: "
      << accumulate(count.begin(), count.end(), 0u, plus<unsigned int>())
    );
}

}} // namespace algorithm::host

} // namespace halmd

#endif /* ! HALMD_ALGORITHM_HOST_PICK_LATTICE_POINTS_HPP */
