/*
 * Copyright © 2012 Peter Colberg
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

#ifndef HALMD_ALGORITHM_GPU_IOTA_HPP
#define HALMD_ALGORITHM_GPU_IOTA_HPP

#include <halmd/algorithm/gpu/iota_kernel.hpp>

#include <iterator>
#include <type_traits>

namespace halmd {

/**
 * Fill range with sequentially increasing values.
 */
template <typename Iterator>
inline typename std::enable_if<
    std::is_same<
        typename std::iterator_traits<Iterator>::value_type
      , unsigned int
    >::value
    && std::is_convertible<
        typename std::iterator_traits<Iterator>::iterator_category
      , cuda::device_random_access_iterator_tag
    >::value
  , void>::type
iota(
    Iterator const& first
  , Iterator const& last
  , unsigned int value
)
{
    unsigned int constexpr threads = 128;
    unsigned int const size = last - first;
    detail::iota_kernel::get().iota.configure((size + threads - 1) / threads,
        threads);
    detail::iota_kernel::get().iota(&*first, size, value);
    cuda::thread::synchronize();
}

} // namespace halmd

#endif /* ! HALMD_ALGORITHM_GPU_IOTA_HPP */
