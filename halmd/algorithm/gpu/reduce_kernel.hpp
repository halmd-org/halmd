/*
 * Copyright © 2008-2012 Peter Colberg
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

#ifndef HALMD_ALGORITHM_GPU_REDUCE_KERNEL_HPP
#define HALMD_ALGORITHM_GPU_REDUCE_KERNEL_HPP

#include <halmd/utility/gpu/shared_memory.hpp>

#include <cuda_wrapper/cuda_wrapper.hpp>

#include <iterator>
#include <stdexcept>

namespace halmd {
namespace detail {

template <typename accumulator_type>
struct reduction_function
{
    typedef typename accumulator_type::iterator iterator;
    typedef typename std::iterator_traits<iterator>::difference_type difference_type;
    typedef void (type)(iterator const, difference_type, accumulator_type*, accumulator_type);
};

} // namespace detail

template <typename accumulator_type>
struct reduction_kernel
{
    typedef typename detail::reduction_function<accumulator_type>::type function_type;
    cuda::function<function_type> reduce;
    static reduction_kernel kernel;
};

} // namespace halmd

#endif /* ! HALMD_ALGORITHM_GPU_REDUCE_KERNEL_HPP */
