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
template <typename tag, typename output_type>
class reduce : private tag
{
public:
    enum { BLOCKS = gpu::reduce::BLOCKS, THREADS = gpu::reduce::THREADS };

    reduce() : g_block_sum(BLOCKS), h_block_sum(BLOCKS) {}

    /**
     * stream parallel reduction kernel
     */
    template <typename input_type>
    void operator()(cuda::vector<input_type> const& g_in, cuda::stream& stream)
    {
	cuda::configure(BLOCKS, THREADS, stream);
	tag::reduce(g_in, g_block_sum, g_in.size());
	cuda::copy(g_block_sum, h_block_sum, stream);
    }

    /**
     * returns sum after CUDA stream synchronisation
     */
    output_type value()
    {
	return tag::value(h_block_sum);
    }

private:
    cuda::vector<output_type> g_block_sum;
    cuda::host::vector<output_type> h_block_sum;
};

namespace tag {

/**
 * sum
 */
struct sum
{
    template <typename T, typename U>
    void reduce(T const& g_in, U& g_block_sum, unsigned int count)
    {
	gpu::reduce::sum(g_in, g_block_sum, count);
    }

    /**
     * returns sum after CUDA stream synchronisation
     */
    template <typename U>
    typename U::value_type value(U const& sum)
    {
	return std::accumulate(sum.begin(), sum.end(), (typename U::value_type) 0);
    }
};

/**
 * absolute maximum
 */
struct max
{
    template <typename T, typename U>
    void reduce(T const& g_in, U& g_block_max, unsigned int count)
    {
	gpu::reduce::max(g_in, g_block_max, count);
    }

    /**
     * returns max after CUDA stream synchronisation
     */
    template <typename U>
    typename U::value_type value(U const& max)
    {
	return *std::max_element(max.begin(), max.end());
    }
};

} // namespace ljgpu::tag

} // namespace ljgpu

#endif /* ! LJGPU_ALGORITHM_REDUCE_HPP */
