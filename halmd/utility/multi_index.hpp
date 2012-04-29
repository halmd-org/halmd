/*
 * Copyright © 2012  Peter Colberg
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

#ifndef HALMD_UTILITY_MULTI_INDEX_HPP
#define HALMD_UTILITY_MULTI_INDEX_HPP

#include <boost/utility/enable_if.hpp>
#ifndef __CUDACC__
# include <cassert>
#endif

#include <halmd/config.hpp>

namespace halmd {

template <unsigned int dim, typename index_type>
inline HALMD_GPU_ENABLED
typename boost::enable_if_c<(dim == index_type::static_size - 1), typename index_type::value_type>::type
multi_index_to_offset(index_type const& index, index_type const& dims)
{
#ifndef __CUDACC__
    assert(index[dim] < dims[dim]);
#endif
    return index[dim];
}

template <unsigned int dim, typename index_type>
inline HALMD_GPU_ENABLED
typename boost::enable_if_c<(dim < index_type::static_size - 1), typename index_type::value_type>::type
multi_index_to_offset(index_type const& index, index_type const& dims)
{
#ifndef __CUDACC__
    assert(index[dim] < dims[dim]);
#endif
    return index[dim] + dims[dim] * multi_index_to_offset<dim + 1>(index, dims);
}

/**
 * Convert multi-dimensional index to one-dimensional offset.
 */
template <typename index_type>
inline HALMD_GPU_ENABLED typename index_type::value_type
multi_index_to_offset(index_type const& index, index_type const& dims)
{
    return multi_index_to_offset<0>(index, dims);
}

template <unsigned int dim, typename index_type>
inline HALMD_GPU_ENABLED
typename boost::enable_if_c<(dim == index_type::static_size - 1), void>::type
offset_to_multi_index(index_type& index, typename index_type::value_type const& offset, index_type const& dims)
{
    index[dim] = offset;
#ifndef __CUDACC__
    assert(index[dim] < dims[dim]);
#endif
}

template <unsigned int dim, typename index_type>
inline HALMD_GPU_ENABLED
typename boost::enable_if_c<(dim < index_type::static_size - 1), void>::type
offset_to_multi_index(index_type& index, typename index_type::value_type const& offset, index_type const& dims)
{
    index[dim] = offset % dims[dim];
#ifndef __CUDACC__
    assert(index[dim] < dims[dim]);
#endif
    offset_to_multi_index<dim + 1>(index, offset / dims[dim], dims);
}

/**
 * Convert one-dimensional offset to multi-dimensional index.
 */
template <typename index_type>
inline HALMD_GPU_ENABLED index_type
offset_to_multi_index(typename index_type::value_type const& offset, index_type const& dims)
{
    index_type index;
    offset_to_multi_index<0>(index, offset, dims);
    return index;
}

} // namespace halmd

#endif /* ! HALMD_UTILITY_MULTI_INDEX_HPP */
