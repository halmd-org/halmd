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

#ifndef HALMD_OBSERVABLES_UTILITY_WAVEVECTORS_HPP
#define HALMD_OBSERVABLES_UTILITY_WAVEVECTORS_HPP

#include <algorithm>
#include <boost/foreach.hpp>
#include <cmath>
#include <map>
#include <numeric>
#include <vector>

#include <halmd/io/logger.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/numeric/gcd.hpp>

namespace halmd { namespace observables { namespace utility
{

/**
 *  construct shells of wavevectors from reciprocal lattice
 *
 * @param wavenumbers:  list of wavenumbers, defines radii of shells
 * @param tolerance:    relative tolerance on wavenumbers, defines thickness of shells
 * @param miller_max:   maximal sum of Miller indices, limiting the density of reciprocal space
 * @returns list of wavevector shells, the wavevectors of each shell have the same magnitude
 *          (within the given tolerance)
 *
 * The constructed wavevectors are of the form @f$ \vec q = (2\pi/L) n (h, k, l) @f$
 * where @f$ 0 \leq h+k+l \leq \text{miller_max} @f$ are coprime and @f$ n @f$ is
 * integer.
 */
template <typename float_type, typename vector_type>
std::multimap<float_type, vector_type> construct_wavevector_shells(
    std::vector<float_type> const& wavenumbers
  , vector_type const& box_length
  , float_type rel_error
  , unsigned int miller_max
)
{
    using namespace std;

    enum { dimension = vector_type::static_size };

    typedef fixed_vector<unsigned int, dimension> index_type;
    typedef multimap<float_type, vector_type> map_type;
    map_type wavevectors;

    LOG_DEBUG("generate wavevectors with a maximum sum of Miller indices: " << miller_max);

    // set of orthogonal basis vectors (stored in a single vector, 2π / L_i)
    vector_type q_basis = element_div(vector_type(2 * M_PI), box_length);

    // generate all wavevectors bounded by q_max and miller_max
    // loop over tuple index (hkl), or (hk) in two dimensions
    index_type hkl(0);
    bool carry_flag = false;
    while (!carry_flag) {
        if (// 1st condition: h+k+l ≤ miller_max
            accumulate(hkl.begin(), hkl.end(), 0u, plus<unsigned int>()) <= miller_max
            // 2nd condition: tuple elements are coprime
         && is_coprime(hkl))
        {
            // construct and store wavevectors along the direction (hkl)
            // with magnitude close to the desired values
            vector_type q0 = element_prod(q_basis, static_cast<vector_type>(hkl));
            float_type q0_norm = norm_2(q0);
            BOOST_FOREACH (float_type q, wavenumbers) {
                // find integer n such that abs(norm_2(n * q0) - q) / q < q_error
                // 1) round to nearest integer
                unsigned int n = floor(q / q0_norm + float_type(.5));
                // 2) check if this is good enough
                if (n > 0 && abs(n * q0_norm / q - 1) < rel_error) {
                    wavevectors.insert(make_pair(q, n * q0));
                }
            }
        }
        // increment tuple (hkl) at end of loop
        carry_flag = true;                   //< assume carry over for 'digit' #0
        for (unsigned int j = 0; carry_flag && j < dimension; ++j) {
            if (hkl[j] < miller_max) {       // miller_max is included in the range of values
                ++hkl[j];                    // increment 'digit'
                carry_flag = false;          //< no carry over, exit loop
            }
            else {
                hkl[j] = 0;                  //< reset this 'digit' and increment next entry
            }
        }
        // carry_flag == true only if carry-over from last 'digit' of hkl
    };
    LOG_TRACE("wavevectors constructed: " << wavevectors.size());

    // output list of wavevectors for debug trace
#ifndef NDEBUG
    BOOST_FOREACH (float_type q, wavenumbers) {
        LOG_TRACE(wavevectors.count(q) << " wavevectors with |q| ≈ " << q);
        typedef pair<typename map_type::const_iterator, typename map_type::const_iterator> range_type;
        for (range_type shell = wavevectors.equal_range(q); shell.first != shell.second; ++shell.first) {
            vector_type const& q_vector = shell.first->second;
            index_type hkl = static_cast<index_type>(round(element_div(q_vector, q_basis)));
            hkl /= greatest_common_divisor(hkl);
            LOG_TRACE("  |q| = " << norm_2(q_vector) << ", (hkl) = " << hkl);
        }
    }
#endif

    return wavevectors;
}

}}} // namespace halmd::observables::utility

#endif /* ! HALMD_OBSERVABLES_UTILITY_WAVEVECTORS_HPP */
