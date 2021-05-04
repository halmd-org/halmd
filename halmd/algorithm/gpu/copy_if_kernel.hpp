/*
 * Copyright © 2015 Nicolas Höft
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

#ifndef HALMD_ALGORITHM_GPU_COPY_IF_KERNEL_HPP
#define HALMD_ALGORITHM_GPU_COPY_IF_KERNEL_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <boost/version.hpp>

#if defined(__CUDACC__) && BOOST_VERSION >= 106900
# define BOOST_NO_CXX11_CONSTEXPR
#endif
#include <boost/function.hpp>
#if defined(__CUDACC__) && BOOST_VERSION >= 106900
# undef BOOST_NO_CXX11_CONSTEXPR
#endif

namespace halmd {
namespace algorithm {
namespace gpu {

/**
 * CUDA C++ wrapper
 */
template<typename InputIterator, typename OutputIterator, typename Predicate>
struct copy_if_wrapper
{
    boost::function<unsigned int (
        InputIterator          // input array
      , unsigned int           // array length
      , Predicate
      , OutputIterator         // output array
    )> copy_if;
    static copy_if_wrapper const kernel;
};

} // namespace gpu
} // namespace algorithm
} // namespace halmd

#endif /* ! HALMD_ALGORITHM_GPU_COPY_IF_KERNEL_HPP */
