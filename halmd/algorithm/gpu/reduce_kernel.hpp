/*
 * Copyright © 2008-2012 Peter Colberg
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

template <unsigned int threads, typename accumulator_type>
struct reduction_kernel
{
    typedef typename reduction_function<accumulator_type>::type function_type;
    cuda::function<function_type> reduce;
    static reduction_kernel const kernel;
};

} // namespace detail

template <typename accumulator_type, unsigned int max_threads = shared_memory_max_threads<accumulator_type>::value>
struct reduction_kernel;

template <typename accumulator_type>
struct reduction_kernel<accumulator_type, 1024>
{
    typedef typename detail::reduction_function<accumulator_type>::type function_type;

    static cuda::function<function_type> const& reduce(unsigned int threads)
    {
        switch (threads) {
          case 1024:
            return detail::reduction_kernel<1024, accumulator_type>::kernel.reduce;
          case 512:
            return detail::reduction_kernel<512, accumulator_type>::kernel.reduce;
          case 256:
            return detail::reduction_kernel<256, accumulator_type>::kernel.reduce;
          case 128:
            return detail::reduction_kernel<128, accumulator_type>::kernel.reduce;
          case 64:
            return detail::reduction_kernel<64, accumulator_type>::kernel.reduce;
          case 32:
            return detail::reduction_kernel<32, accumulator_type>::kernel.reduce;
          default:
            throw std::logic_error("invalid number of reduction threads");
        }
    }
};

template <typename accumulator_type>
struct reduction_kernel<accumulator_type, 512>
{
    typedef typename detail::reduction_function<accumulator_type>::type function_type;

    static cuda::function<function_type> const& reduce(unsigned int threads)
    {
        switch (threads) {
          case 512:
            return detail::reduction_kernel<512, accumulator_type>::kernel.reduce;
          case 256:
            return detail::reduction_kernel<256, accumulator_type>::kernel.reduce;
          case 128:
            return detail::reduction_kernel<128, accumulator_type>::kernel.reduce;
          case 64:
            return detail::reduction_kernel<64, accumulator_type>::kernel.reduce;
          case 32:
            return detail::reduction_kernel<32, accumulator_type>::kernel.reduce;
          default:
            throw std::logic_error("invalid number of reduction threads");
        }
    }
};

template <typename accumulator_type>
struct reduction_kernel<accumulator_type, 256>
{
    typedef typename detail::reduction_function<accumulator_type>::type function_type;

    static cuda::function<function_type> const& reduce(unsigned int threads)
    {
        switch (threads) {
          case 256:
            return detail::reduction_kernel<256, accumulator_type>::kernel.reduce;
          case 128:
            return detail::reduction_kernel<128, accumulator_type>::kernel.reduce;
          case 64:
            return detail::reduction_kernel<64, accumulator_type>::kernel.reduce;
          case 32:
            return detail::reduction_kernel<32, accumulator_type>::kernel.reduce;
          default:
            throw std::logic_error("invalid number of reduction threads");
        }
    }
};

template <typename accumulator_type>
struct reduction_kernel<accumulator_type, 128>
{
    typedef typename detail::reduction_function<accumulator_type>::type function_type;

    static cuda::function<function_type> const& reduce(unsigned int threads)
    {
        switch (threads) {
          case 128:
            return detail::reduction_kernel<128, accumulator_type>::kernel.reduce;
          case 64:
            return detail::reduction_kernel<64, accumulator_type>::kernel.reduce;
          case 32:
            return detail::reduction_kernel<32, accumulator_type>::kernel.reduce;
          default:
            throw std::logic_error("invalid number of reduction threads");
        }
    }
};

template <typename accumulator_type>
struct reduction_kernel<accumulator_type, 64>
{
    typedef typename detail::reduction_function<accumulator_type>::type function_type;

    static cuda::function<function_type> const& reduce(unsigned int threads)
    {
        switch (threads) {
          case 64:
            return detail::reduction_kernel<64, accumulator_type>::kernel.reduce;
          case 32:
            return detail::reduction_kernel<32, accumulator_type>::kernel.reduce;
          default:
            throw std::logic_error("invalid number of reduction threads");
        }
    }
};

template <typename accumulator_type>
struct reduction_kernel<accumulator_type, 32>
{
    typedef typename detail::reduction_function<accumulator_type>::type function_type;

    static cuda::function<function_type> const& reduce(unsigned int threads)
    {
        switch (threads) {
          case 32:
            return detail::reduction_kernel<32, accumulator_type>::kernel.reduce;
          default:
            throw std::logic_error("invalid number of reduction threads");
        }
    }
};

} // namespace halmd

#endif /* ! HALMD_ALGORITHM_GPU_REDUCE_KERNEL_HPP */
