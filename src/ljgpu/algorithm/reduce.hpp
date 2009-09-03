/* Parallel reduction kernel
 *
 * Copyright © 2008-2009  Peter Colberg
 *
 * This program is free software: you can redistribute it and/or modify
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

#ifndef LJGPU_ALGORITHM_REDUCE_HPP
#define LJGPU_ALGORITHM_REDUCE_HPP

#include <algorithm>
#include <cuda_wrapper.hpp>
#include <ljgpu/algorithm/gpu/reduce.hpp>
#include <numeric>

namespace ljgpu
{

/*
 * Parallel reduction
 */
template <
    template <typename> class tag,
    typename gpu_output_type,
    typename output_type = gpu_output_type>
class reduce
{
public:
    enum { BLOCKS = gpu::reduce::BLOCKS, THREADS = gpu::reduce::THREADS };

    reduce() : g_block_sum(BLOCKS), h_block_sum(BLOCKS) {}

    /**
     * stream parallel reduction kernel
     */
    template <typename gpu_input_type>
    void operator()(cuda::vector<gpu_input_type> const& g_in, cuda::stream& stream)
    {
        cuda::configure(BLOCKS, THREADS, stream);
        tag<output_type>::reduce(g_in, g_block_sum, g_in.size());
        cuda::copy(g_block_sum, h_block_sum, stream);
    }

    /**
     * returns sum after CUDA stream synchronisation
     */
    output_type value()
    {
        return tag<output_type>::value(h_block_sum);
    }

private:
    cuda::vector<gpu_output_type> g_block_sum;
    cuda::host::vector<gpu_output_type> h_block_sum;
};

namespace tag {

/**
 * sum
 */
template <typename output_type>
struct sum
{
    template <typename gpu_input_type, typename gpu_output_type>
    static void reduce(cuda::vector<gpu_input_type> const& g_in,
                       cuda::vector<gpu_output_type>& g_block_sum, unsigned int count)
    {
        gpu::reduce::sum(g_in, g_block_sum, count);
    }

    /**
     * returns sum after CUDA stream synchronisation
     */
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
    static void reduce(cuda::vector<gpu_input_type> const& g_in,
                       cuda::vector<gpu_output_type>& g_block_sum, unsigned int count)
    {
        gpu::reduce::sum_of_squares(g_in, g_block_sum, count);
    }

    /**
     * returns sum after CUDA stream synchronisation
     */
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
    static void reduce(cuda::vector<gpu_input_type> const& g_in,
                       cuda::vector<gpu_output_type>& g_block_max, unsigned int count)
    {
        gpu::reduce::max(g_in, g_block_max, count);
    }

    /**
     * returns max after CUDA stream synchronisation
     */
    template <typename gpu_output_type>
    static output_type value(cuda::host::vector<gpu_output_type> const& max)
    {
        return *std::max_element(max.begin(), max.end());
    }
};

} // namespace ljgpu::tag

} // namespace ljgpu

#endif /* ! LJGPU_ALGORITHM_REDUCE_HPP */
