/*
 * Copyright © 2008-2012  Peter Colberg
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

#ifndef HALMD_ALGORITHM_GPU_REDUCE_HPP
#define HALMD_ALGORITHM_GPU_REDUCE_HPP

#include <algorithm> // std::for_each
#include <boost/utility.hpp> // boost::noncopyable
#include <cassert> // assert
#include <stdexcept>

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <halmd/algorithm/gpu/reduce_kernel.hpp>

namespace halmd {

/**
 * Reduce an input array on GPU using an accumulator.
 */
template <
    typename accumulator_type
  , unsigned int max_threads = shared_memory_max_threads<accumulator_type>::value
  , unsigned int arity = detail::reduction_arity<accumulator_type>::value
>
class reduction;

template <typename accumulator_type, unsigned int max_threads>
class reduction<accumulator_type, max_threads, 1>
  : boost::noncopyable
{
private:
    typedef reduction_kernel<accumulator_type, max_threads> kernel_type;
    typedef typename kernel_type::function_type function_type;
    typedef typename accumulator_type::argument_type argument_type;

public:
    /**
     * Allocate reduction buffers in GPU and host memory.
     *
     * @param blocks number of blocks in execution grid
     * @param threads number of threads per block
     *
     * The number of threads per block must be power of 2, and at least 32
     * (size of a wrap). CUDA PTX < 2.0 allows a maximum number of threads
     * of 512, while CUDA PTX ≥ 2.0 allows 1024 threads per block.
     */
    reduction(unsigned int blocks = 16, unsigned int threads = max_threads);

    /**
     * Reduce values in input array using given accumulator.
     *
     * @param g_input input array in GPU memory
     * @param acc reduction accumulator
     */
    template <typename InputIterator>
    accumulator_type operator()(
        InputIterator first
      , InputIterator last
      , accumulator_type const& acc = accumulator_type()
    );

private:
    /** kernel execution parameters */
    cuda::config dim_;
    /** accumulators per block in GPU memory */
    cuda::vector<accumulator_type> g_block_;
    /** accumulators per block in pinned host memory */
    cuda::host::vector<accumulator_type> h_block_;
    /** reduce kernel function matching number of threads per block */
    cuda::function<function_type> reduce_;
};

template <typename accumulator_type, unsigned int max_threads>
inline reduction<accumulator_type, max_threads, 1>::reduction(
    unsigned int blocks
  , unsigned int threads
)
  : dim_(blocks, threads)
  , g_block_(blocks)
  , reduce_(kernel_type::reduce(threads))
{
    // avoid DefaultConstructible requirement on accumulator_type
    h_block_.reserve(blocks);
}

template <typename accumulator_type, unsigned int max_threads>
template <typename InputIterator>
inline accumulator_type reduction<accumulator_type, max_threads, 1>::operator()(
    InputIterator first
  , InputIterator last
  , accumulator_type const& acc
)
{
    cuda::configure(dim_.grid, dim_.block);
    reduce_(&*first, last - first, g_block_, acc);
    assert(g_block_.size() == h_block_.capacity());
    cuda::copy(g_block_, h_block_, g_block_.size());
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
template <typename InputIterator, typename accumulator_type>
inline accumulator_type reduce(
    InputIterator first
  , InputIterator last
  , accumulator_type const& acc
)
{
    return reduction<accumulator_type>()(first, last, acc);
}

template <typename accumulator_type, unsigned int max_threads>
class reduction<accumulator_type, max_threads, 2>
  : boost::noncopyable
{
private:
    typedef reduction_kernel<accumulator_type, max_threads> kernel_type;
    typedef typename kernel_type::function_type function_type;
    typedef typename accumulator_type::first_argument_type first_argument_type;
    typedef typename accumulator_type::second_argument_type second_argument_type;

public:
    /**
     * Allocate reduction buffers in GPU and host memory.
     *
     * @param blocks number of blocks in execution grid
     * @param threads number of threads per block
     *
     * The number of threads per block must be power of 2, and at least 32
     * (size of a wrap). CUDA PTX < 2.0 allows a maximum number of threads
     * of 512, while CUDA PTX ≥ 2.0 allows 1024 threads per block.
     */
    reduction(unsigned int blocks = 16, unsigned int threads = max_threads);

    /**
     * Reduce values in input array using binary accumulator.
     *
     * @param g_input input array in GPU memory
     * @param acc reduction accumulator
     */
    template <typename InputIterator1, typename InputIterator2>
    accumulator_type operator()(
        InputIterator1 first1
      , InputIterator1 last1
      , InputIterator2 first2
      , accumulator_type const& acc = accumulator_type()
    );

private:
    /** kernel execution parameters */
    cuda::config dim_;
    /** accumulators per block in GPU memory */
    cuda::vector<accumulator_type> g_block_;
    /** accumulators per block in pinned host memory */
    cuda::host::vector<accumulator_type> h_block_;
    /** reduce kernel function matching number of threads per block */
    cuda::function<function_type> reduce_;
};

template <typename accumulator_type, unsigned int max_threads>
inline reduction<accumulator_type, max_threads, 2>::reduction(
    unsigned int blocks
  , unsigned int threads
)
  : dim_(blocks, threads)
  , g_block_(blocks)
  , reduce_(kernel_type::reduce(threads))
{
    // avoid DefaultConstructible requirement on accumulator_type
    h_block_.reserve(blocks);
}

template <typename accumulator_type, unsigned int max_threads>
template <typename InputIterator1, typename InputIterator2>
inline accumulator_type reduction<accumulator_type, max_threads, 2>::operator()(
    InputIterator1 first1
  , InputIterator1 last1
  , InputIterator2 first2
  , accumulator_type const& acc
)
{
    cuda::configure(dim_.grid, dim_.block);
    reduce_(&*first1, &*first2, last1 - first1, g_block_, acc);
    assert(g_block_.size() == h_block_.capacity());
    cuda::copy(g_block_, h_block_, g_block_.size());
    return std::for_each(h_block_.begin(), h_block_.begin() + h_block_.capacity(), acc);
}

/**
 * Reduce values of two input arrays using binary accumulator.
 *
 * @param g_first first input array in GPU memory
 * @param g_second second input array in GPU memory
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
template <typename InputIterator1, typename InputIterator2, typename accumulator_type>
inline accumulator_type reduce(
    InputIterator1 first1
  , InputIterator1 last1
  , InputIterator2 first2
  , accumulator_type const& acc
)
{
    return reduction<accumulator_type>()(first1, last1, first2, acc);
}

} // namespace halmd

#endif /* ! HALMD_ALGORITHM_GPU_REDUCE_HPP */
