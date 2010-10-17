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

#ifndef HALMD_ALGORITHM_HOST_PERMUTE_HPP
#define HALMD_ALGORITHM_HOST_PERMUTE_HPP

#include <vector>

namespace halmd
{
namespace algorithm { namespace host
{

/*
 * Permute sequence in-place according to integer index sequence
 *
 * This is a faster variant of the algorithm described in
 *
 *   Donald E. Knuth, Selected Papers on Analysis of Algorithms
 *   (Stanford, California: Center for the Study of Language and
 *   Information - CSLI Lecture Notes, no. 102), 2000.
 */
template <typename input_iterator, typename index_iterator>
void permute(input_iterator first, input_iterator last, index_iterator index)
{
    typedef typename input_iterator::difference_type difference_type;
    typedef typename input_iterator::value_type value_type;

    std::vector<bool> follower(last - first, false);
    for (difference_type i = 0; i < last - first; ++i) {
        if (!follower[i]) {
            value_type temp = *(first + i);
            difference_type j = i;
            for (difference_type k = *(index + j); k != i; j = k, k = *(index + j)) {
                *(first + j) = *(first + k);
                follower[k] = true;
            }
            *(first + j) = temp;
        }
    }
}

}} // namespace algorithm::host

} // namespace halmd

#endif /* ! HALMD_ALGORITHM_HOST_PERMUTE_HPP */
