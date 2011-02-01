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
 * @param tolerance:    relative tolerance on wavenumbers, defines thickness of shells
 * @param miller_max:   maximal Miller index, limiting the density of reciprocal space
 * @returns list of wavevector shells, the wavevectors of each shell have the same magnitude
 *          (within the given tolerance)
 *
 * The constructed wavevectors are of the form @f$ \vec q = (2\pi/L) n (h, k, l) @f$
 * where @f$ 0 \leq h,k,l \leq \text{miller_max} @f$ are coprime and @f$ n @f$ is
 * integer.
 */
template <typename float_type, typename vector_type>
std::vector<std::vector<vector_type> >
construct_wavevector_shells(
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
    map_type q_list;
    vector<vector<vector_type> > q_shells;

    // maximum wavenumber including tolerance
    float_type q_max = *max_element(wavenumbers.begin(), wavenumbers.end()) * (1 + rel_error);
    float_type qq_max = q_max * q_max;

    // set of orthogonal basis vectors (stored in a single vector, 2π / L_i)
    vector_type q_basis = element_div(vector_type(2 * M_PI), box_length);

    // reduce maximum Miller index if q_max is not sufficiently large
    // *min_element(q_basis) * miller_max ≤ q_max
    miller_max = min(
        miller_max
      , static_cast<unsigned int>(floor(q_max / *min_element(q_basis.begin(), q_basis.end())))
    );
    LOG_DEBUG("generate wavevectors with maximum Miller index: " << miller_max);

    // generate all wavevectors bounded by q_max and miller_max
    // loop over tuple index (hkl), or (hk) in two dimensions
    index_type hkl(0);
    bool carry_flag = false;
    while (!carry_flag) {
        // choose coprime tuples only
        if (is_coprime(hkl)) {
            // construct and store wavevectors along the direction (hkl)
            // with magnitude less or equal to q_max
            vector_type q = element_prod(q_basis, static_cast<vector_type>(hkl));
            float_type qq = inner_prod(q, q);
            for (unsigned int n = 1; n * n * qq < qq_max; ++n) {
                q_list.insert(make_pair(n * n * qq, n * q));
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
    LOG_TRACE("wavevectors constructed: " << q_list.size());

    // sort wavevectors into shells
    q_shells.resize(wavenumbers.size());
    for (unsigned int i = 0; i < wavenumbers.size(); ++i) {
        // find wavevectors with magnitude close to given wavenumber
        typename map_type::iterator first, last, it;
        first = q_list.lower_bound(pow(wavenumbers[i] * (1 - rel_error), 2));
        last  = q_list.upper_bound(pow(wavenumbers[i] * (1 + rel_error), 2));
        for (it = first; it != last; ++it) {
            q_shells[i].push_back(it->second);
        }
        q_list.erase(first, last);
    }

    // output list of wavevectors for debug trace
#ifndef NDEBUG
    LOG_TRACE("wavevectors dropped: " << q_list.size());
    for (unsigned int i = 0; i < q_shells.size(); ++i) {
        vector<vector_type> const& shell = q_shells[i];
        LOG_TRACE(shell.size() << " wavevectors with |q| ≈ " << wavenumbers[i]);
        for (unsigned int j = 0; j < shell.size(); ++j) {
            vector_type const& q = shell[j];
            index_type hkl = static_cast<index_type>(round(element_div(q, q_basis)));
            hkl /= greatest_common_divisor(hkl);
            LOG_TRACE("  |q| = " << norm_2(q) << ", (hkl) = " << hkl);
        }
    }
#endif

    return q_shells;
}

}}} // namespace halmd::observables::utility

#endif /* ! HALMD_OBSERVABLES_UTILITY_WAVEVECTORS_HPP */
