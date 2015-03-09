/*
 * Copyright © 2015 Nicolas Höft
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

#ifndef HALMD_ALGORITHM_GPU_COPY_IF_HPP
#define HALMD_ALGORITHM_GPU_COPY_IF_HPP

#include <iterator>
#include <type_traits>

#include <halmd/algorithm/gpu/copy_if_kernel.hpp>

namespace halmd {

/**
 * Copy elements in the range [first, last) to output, if the predicate is
 * fulfilled.
 * Predicate is a functor with signature bool operator()(T const& value)
 *
 * Note that this function is designed to work with cuda::vector both as input
 * and output iterator
 */
template <typename Iterator, typename Predicate>
inline typename std::enable_if<
    std::is_convertible<
        typename std::iterator_traits<Iterator>::iterator_category
      , cuda::device_random_access_iterator_tag
    >::value
  , Iterator>::type
copy_if(
    Iterator const& first
  , Iterator const& last
  , Iterator const& output
  , Predicate pred
)
{
    typedef typename std::iterator_traits<Iterator>::pointer pointer;

    // the iterators from cuda::vector are not supported in CUDA kernels,
    // therefore use pointers of cuda::vector elements instead of iterators
    auto const& kernel = halmd::algorithm::gpu::copy_if_wrapper<pointer, pointer, Predicate>::kernel;
    unsigned int output_size = kernel.copy_if(
        &*first
      , last - first
      , pred
      , &*output
    );

    return output + output_size;
}

} // namespace halmd

#endif /* ! HALMD_ALGORITHM_GPU_COPY_IF_HPP */
