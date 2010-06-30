/*
 * Copyright © 2008-2009  Peter Colberg
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

#include <algorithm>
#include <numeric>

#include <cuda_wrapper.hpp>
#include <halmd/algorithm/gpu/reduce_kernel.cuh>

namespace halmd
{
namespace algorithm { namespace gpu
{

/*
 * Parallel reduction
 */
template <
    template <typename> class tag
  , typename gpu_output_type
  , typename output_type = gpu_output_type
  , int blocks = 16
  , int threads = (64 << DEVICE_SCALE)
>
struct reduce
{
    template <typename gpu_input_type>
    output_type operator()(cuda::vector<gpu_input_type> const& g_in)
    {
        cuda::vector<gpu_output_type> g_block_sum(blocks);
        cuda::host::vector<gpu_output_type> h_block_sum(blocks);
        cuda::configure(blocks, threads);
        tag<output_type>::reduce(g_in, g_block_sum, g_in.size());
        cuda::copy(g_block_sum, h_block_sum);
        return tag<output_type>::value(h_block_sum);
    }
};

namespace tag
{

/**
 * sum
 */
template <typename output_type>
struct sum
{
    template <typename gpu_input_type, typename gpu_output_type>
    static void reduce(
        cuda::vector<gpu_input_type> const& g_in
      , cuda::vector<gpu_output_type>& g_block_sum
      , unsigned int count
    )
    {
        reduce_wrapper<threads, gpu_input_type, gpu_output_type>::sum(g_in, g_block_sum, count);
    }

    template <typename gpu_output_type>
    static output_type value(cuda::host::vector<gpu_output_type> const& sum)
    {
        return std::accumulate(sum.begin(), sum.end(), output_type(0));
    }
};

/**
 * sum of squares
 */
template <typename output_type>
struct sum_of_squares
{
    template <typename gpu_input_type, typename gpu_output_type>
    static void reduce(
        cuda::vector<gpu_input_type> const& g_in
      , cuda::vector<gpu_output_type>& g_block_sum
      , unsigned int count
    )
    {
        reduce_wrapper<threads, gpu_input_type, gpu_output_type>::sum_of_squares(g_in, g_block_sum, count);
    }

    template <typename gpu_output_type>
    static output_type value(cuda::host::vector<gpu_output_type> const& sum)
    {
        return std::accumulate(sum.begin(), sum.end(), output_type(0));
    }
};

/**
 * absolute maximum
 */
template <typename output_type>
struct max
{
    template <typename gpu_input_type, typename gpu_output_type>
    static void reduce(
        cuda::vector<gpu_input_type> const& g_in
      , cuda::vector<gpu_output_type>& g_block_max
      , unsigned int count
    )
    {
        reduce_wrapper<threads, gpu_input_type, gpu_output_type>::max(g_in, g_block_max, count);
    }

    template <typename gpu_output_type>
    static output_type value(cuda::host::vector<gpu_output_type> const& max)
    {
        return *std::max_element(max.begin(), max.end());
    }
};

} // namespace tag

}} // namespace algorithm::gpu

} // namespace halmd

#endif /* ! HALMD_ALGORITHM_GPU_REDUCE_HPP */
