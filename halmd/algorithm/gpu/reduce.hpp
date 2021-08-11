/*
 * Copyright © 2008-2012 Peter Colberg
 * Copyright © 2021      Jaslo Ziska
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

#ifndef HALMD_ALGORITHM_GPU_REDUCE_HPP
#define HALMD_ALGORITHM_GPU_REDUCE_HPP

#include <halmd/algorithm/gpu/reduce_kernel.hpp>
#include <halmd/utility/gpu/configure_kernel.hpp>

#include <cuda_wrapper/cuda_wrapper.hpp>

#include <algorithm>
#include <cassert>

namespace halmd {

/**
 * Reduce an input array on GPU using an accumulator.
 */
template <typename accumulator_type>
class reduction
{
private:
    typedef reduction_kernel<accumulator_type> kernel_type;

public:
    /**
     * Allocate reduction buffers in GPU and host memory.
     *
     * @param blocks number of blocks in execution grid
     * @param threads number of threads per block
     *
     * The number of threads per block must be a multiple of 32
     * (size of a warp). CUDA PTX ≥ 2.0 allows up to 1024 threads per block.
     */
    reduction(unsigned int blocks = 16, unsigned int threads = 1024);

    /**
     * Reduce values in input array using given accumulator.
     *
     * @param g_input input array in GPU memory
     * @param acc reduction accumulator
     */
    accumulator_type operator()(
        typename accumulator_type::iterator first
      , typename accumulator_type::iterator last
      , accumulator_type const& acc = accumulator_type()
    );

    /** deleted implicit copy constructor */
    reduction(reduction const&) = delete;
    /** deleted implicit assignment operator */
    reduction& operator=(reduction const&) = delete;

private:
    /** accumulators per block in GPU memory */
    cuda::memory::device::vector<accumulator_type> g_block_;
    /** accumulators per block in pinned host memory */
    cuda::memory::host::vector<accumulator_type> h_block_;
};

template <typename accumulator_type>
inline reduction<accumulator_type>::reduction(
    unsigned int blocks
  , unsigned int threads
)
{
    cuda::config dim = configure_kernel(kernel_type::kernel.reduce, cuda::config(blocks, threads), false);
    g_block_.resize(dim.grid.x);
    h_block_.reserve(dim.grid.x); // avoid DefaultConstructible requirement on accumulator_type
}

template <typename accumulator_type>
inline accumulator_type reduction<accumulator_type>::operator()(
    typename accumulator_type::iterator first
  , typename accumulator_type::iterator last
  , accumulator_type const& acc
)
{
    kernel_type::kernel.reduce(first, last - first, g_block_, acc);
    assert(g_block_.size() == h_block_.capacity());
    cuda::copy(g_block_.begin(), g_block_.end(), h_block_.begin());
    return std::for_each(h_block_.begin(), h_block_.begin() + h_block_.capacity(), acc);
}

/**
 * Reduce values in input array using unary accumulator.
 *
 * @param g_input input array in GPU memory
 * @param acc reduction accumulator
 *
 * This function is provided for convenience in unit tests only.
 *
 * You are strongly advised to pre-allocate the reduction buffers by
 * defining a reduction<accumulator_type> member in your C++ class,
 * and repeatedly invoke this functor to reduce values.
 *
 * The allocation of pinned host memory may be unpredictably slow, and take
 * longer than the reduction itself for small arrays. For details, see and
 * run the reduction unit test, and compare the “global” and “local”
 * benchmark results.
 */
template <typename accumulator_type>
inline accumulator_type reduce(
    typename accumulator_type::iterator first
  , typename accumulator_type::iterator last
  , accumulator_type const& acc
)
{
    return reduction<accumulator_type>()(first, last, acc);
}

} // namespace halmd

#endif /* ! HALMD_ALGORITHM_GPU_REDUCE_HPP */
