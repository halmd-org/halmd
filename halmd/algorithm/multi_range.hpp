/*
 * Copyright © 2011  Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#ifndef HALMD_ALGORITHM_MULTI_RANGE_HPP
#define HALMD_ALGORITHM_MULTI_RANGE_HPP

#include <boost/utility/enable_if.hpp>

#include <halmd/config.hpp>

namespace halmd {

template <size_t dimension, typename index_type, typename function_type>
inline HALMD_GPU_ENABLED
typename boost::enable_if_c<dimension == index_type::static_size, void>::type
multi_range_for_each(
    index_type const& it
  , index_type const&
  , index_type const&
  , function_type& f
)
{
    f(it);
}

template <size_t dimension, typename index_type, typename function_type>
inline HALMD_GPU_ENABLED
typename boost::enable_if_c<dimension < index_type::static_size, void>::type
multi_range_for_each(
    index_type& it
  , index_type const& first
  , index_type const& last
  , function_type& f
)
{
    it[dimension] = first[dimension];
    while (it[dimension] != last[dimension]) {
        multi_range_for_each<dimension + 1>(it, first, last, f);
        ++it[dimension];
    }
}

/**
 * Apply function to each index in multi-dimensional range [first, last)
 *
 * @returns index to end of range
 */
template <typename index_type, typename function_type>
inline HALMD_GPU_ENABLED
function_type multi_range_for_each(
    index_type const& first
  , index_type const& last
  , function_type f
)
{
    index_type it(first);
    multi_range_for_each<0>(it, first, last, f);
    return f;
}

template <size_t dimension, typename index_type, typename predicate_type>
inline HALMD_GPU_ENABLED
typename boost::enable_if_c<dimension == index_type::static_size, bool>::type
multi_range_find_if(
    index_type const& it
  , index_type const&
  , index_type const&
  , predicate_type& pred
)
{
    return pred(it);
}

template <size_t dimension, typename index_type, typename predicate_type>
inline HALMD_GPU_ENABLED
typename boost::enable_if_c<dimension < index_type::static_size, bool>::type
multi_range_find_if(
    index_type& it
  , index_type const& first
  , index_type const& last
  , predicate_type& pred
)
{
    it[dimension] = first[dimension];
    while (it[dimension] != last[dimension]) {
        if (multi_range_find_if<dimension + 1>(it, first, last, pred)) {
            return true;
        }
        ++it[dimension];
    }
    return false;
}

/**
 * Find first index in multi-dimensional range for which predicate yields true
 *
 * @returns first matching index, or index to end of range otherwise
 */
template <typename index_type, typename predicate_type>
inline HALMD_GPU_ENABLED
index_type multi_range_find_if(
    index_type const& first
  , index_type const& last
  , predicate_type pred
)
{
    index_type it(first);
    multi_range_find_if<0>(it, first, last, pred);
    return it;
}

} // namespace halmd

#endif /* HALMD_ALGORITHM_MULTI_RANGE_HPP */
