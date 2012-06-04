/*
 * Copyright © 2008-2010, 2012 Peter Colberg
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

#include <cassert>
#include <iterator>
#include <limits>
#include <type_traits>
#include <vector>

namespace halmd {

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
typename std::enable_if<
    std::is_convertible<
        typename std::iterator_traits<input_iterator>::iterator_category
      , std::random_access_iterator_tag
    >::value
    && std::is_convertible<
        typename std::iterator_traits<index_iterator>::iterator_category
      , std::random_access_iterator_tag
    >::value
    && std::numeric_limits<
        typename std::iterator_traits<index_iterator>::value_type
    >::is_integer
  , index_iterator>::type
permute(
    input_iterator const& first
  , input_iterator const& last
  , index_iterator const& index
)
{
    typedef typename std::iterator_traits<input_iterator>::difference_type difference_type;
    typedef typename std::iterator_traits<input_iterator>::value_type value_type;

    difference_type const count = last - first;
    std::vector<bool> follower(count, false);
    for (difference_type i = 0; i < count; ++i) {
        if (!follower[i]) {
            value_type temp = *(first + i);
            difference_type j = i;
            for (difference_type k = *(index + j); k != i; j = k, k = *(index + j)) {
                assert(k < count);
                *(first + j) = *(first + k);
                follower[k] = true;
            }
            *(first + j) = temp;
        }
    }
    return index + count;
}

} // namespace halmd

#endif /* ! HALMD_ALGORITHM_HOST_PERMUTE_HPP */
