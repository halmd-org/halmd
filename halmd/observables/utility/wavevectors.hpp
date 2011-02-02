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
 * @param box_length:   lengths of box edges, defines granularity of reciprocal space
 * @param tolerance:    relative tolerance on wavenumbers, defines thickness of shells
 * @param max_count:    maximal number of wavevectors in each shell
 * @returns list of wavevector shells, the wavevectors of each shell have the same magnitude
 *          (within the given tolerance)
 *
 * The constructed wavevectors are of the form @f$ \vec q = (2\pi/L) n (h, k, l) @f$
 * where the Miller indices @f$ (hkl) @f$ are coprime and @f$ n @f$ is integer. Wavevectors
 * with a smaller sum @f$ h+k+l @f$ are preferred.
 *
 * Caveat: since the number of wavevectors per shell is sharply limited, the returned
 * wavevectors are not completely symmetric in space (some possible combinations
 * @f$ (h,k,l) @f$ may be missing for the largest sum @f$ h+k+l @f$)
 */
template <typename float_type, typename vector_type>
std::multimap<float_type, vector_type> construct_wavevector_shells(
    std::vector<float_type> const& wavenumbers
  , vector_type const& box_length
  , float_type tolerance
  , unsigned int max_count
)
{
    using namespace std;

    enum { dimension = vector_type::static_size };

    typedef fixed_vector<unsigned int, dimension> index_type;
    typedef multimap<float_type, vector_type> map_type;
    map_type wavevectors;

    // determine maximum Miller index that fits in the wavenumber range
    // (2π / L_i) × (h,k,l) ≤ q_max ⇒ (h,k,l) ≤ q_max max(L_i) / 2π
    float_type q_max = *max_element(wavenumbers.begin(), wavenumbers.end()) * (1 + tolerance);
    unsigned int miller_max = static_cast<unsigned int>(floor(
        q_max * *max_element(box_length.begin(), box_length.end()) / (2 * M_PI)
    ));
    LOG_DEBUG("generate wavevectors with a maximum sum of Miller indices: " << miller_max);

    // set of orthogonal basis vectors (stored in a single vector, 2π / L_i)
    vector_type q_basis = element_div(vector_type(2 * M_PI), box_length);

    // generate all wavevectors bounded h+k+l ≤ miller_max,
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
            // construct and store wavevectors along the direction (hkl)
            // with magnitude close to the desired values
            vector_type q0 = element_prod(q_basis, static_cast<vector_type>(hkl));
            float_type q0_norm = norm_2(q0);
            BOOST_FOREACH (float_type q, wavenumbers) {
                if (wavevectors.count(q) < max_count) {
                    // find integer n such that abs(norm_2(n * q0) - q) / q < tolerance
                    // 1) round to nearest integer
                    unsigned int n = floor(q / q0_norm + float_type(.5));
                    // 2) check if this is good enough
                    if (n > 0 && abs(n * q0_norm - q) < q * tolerance) {
                        wavevectors.insert(make_pair(q, n * q0));
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
    LOG_TRACE("wavevectors constructed: " << wavevectors.size());

    // output list of wavevectors for debug trace
    BOOST_FOREACH (float_type q, wavenumbers) {
        unsigned int count = wavevectors.count(q);
        LOG_TRACE(count << " wavevectors with |q| ≈ " << q);
        if (!count) {
            LOG_WARNING("No wavevector compatible with |q| ≈ " << q);
        }
#ifndef NDEBUG
        typedef pair<typename map_type::const_iterator, typename map_type::const_iterator> range_type;
        for (range_type shell = wavevectors.equal_range(q); shell.first != shell.second; ++shell.first) {
            vector_type const& q_vector = shell.first->second;
            index_type hkl = static_cast<index_type>(round(element_div(q_vector, q_basis)));
            hkl /= greatest_common_divisor(hkl);
            LOG_TRACE("  |q| = " << norm_2(q_vector) << ", (hkl) = " << hkl);
        }
#endif
    }

    return wavevectors;
}

}}} // namespace halmd::observables::utility

#endif /* ! HALMD_OBSERVABLES_UTILITY_WAVEVECTORS_HPP */
