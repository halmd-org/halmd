/*
 * Copyright © 2008, 2012 Peter Colberg
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

#ifndef HALMD_ALGORITHM_HOST_RADIX_SORT_HPP
#define HALMD_ALGORITHM_HOST_RADIX_SORT_HPP

#include <algorithm>
#include <array>
#include <iterator>
#include <limits>
#include <vector>

namespace halmd {

/**
 * In-place radix sort.
 *
 * Refer to the unit test for a performance comparison with std::sort.
 */
template <typename Iterator>
typename std::enable_if<
    std::is_convertible<
        typename std::iterator_traits<Iterator>::iterator_category
      , std::forward_iterator_tag
    >::value
    && std::numeric_limits<
        typename std::iterator_traits<Iterator>::value_type
    >::is_integer
  , void>::type radix_sort(Iterator const& first, Iterator const& last)
{
    typedef typename std::iterator_traits<Iterator>::value_type value_type;
    typedef std::vector<value_type> bucket_type;

    unsigned int constexpr digits = std::numeric_limits<value_type>::digits;
    unsigned int constexpr shift = 8;
    unsigned int constexpr size = 1 << shift;
    unsigned int constexpr mask = size - 1;

    std::array<bucket_type, size> buckets;
    for (unsigned int offset = 0; offset < digits; offset += shift) {
        std::for_each(
            first
          , last
          , [&](value_type const& value) {
                buckets[(value >> offset) & mask].push_back(value);
            }
        );
        Iterator output = first;
        for (bucket_type& bucket : buckets) {
            output = std::copy(
                bucket.begin()
              , bucket.end()
              , output
            );
            bucket.clear();
        }
    }
}

} // namespace halmd

#endif /* ! HALMD_ALGORITHM_HOST_RADIX_SORT_HPP */
