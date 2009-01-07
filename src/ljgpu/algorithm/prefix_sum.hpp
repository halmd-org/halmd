/* Parallel exclusive prefix sum
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

#ifndef LJGPU_ALGORITHM_PREFIX_SUM_HPP
#define LJGPU_ALGORITHM_PREFIX_SUM_HPP

#include <ljgpu/algorithm/gpu/prefix_sum.hpp>
#include <stdexcept>
#include <vector>

namespace ljgpu
{

/*
 * Parallel exclusive prefix sum
 */
template <typename T>
class prefix_sum
{
public:
    prefix_sum() : count(0), threads(0) {}

    /**
     * allocate parallel exclusive prefix sum for given element count
     */
    prefix_sum(uint const& count, uint const& threads) : count(count), threads(threads)
    {
	if (threads & (threads - 1)) {
	    throw std::logic_error("prefix sum threads must be a power of 2");
	}

	// compute number of CUDA execution blocks for each level of recursion
	uint n = count;
	do {
	    n = (n + (2 * threads) - 1) / (2 * threads);
	    blocks.push_back(n);
	} while (n > 1);

	// allocate vector with empty GPU global memory vectors
	g_sum.resize(blocks.size());

	// allocate GPU global memory vectors
	for (uint i = 0; i < blocks.size(); ++i) {
	    g_sum[i].resize(blocks[i]);
	}
    }

    /**
     * reallocate parallel exclusive prefix sum for given element count
     */
    void resize(uint const& count, uint const& threads)
    {
	prefix_sum temp(count, threads);
	swap(temp);
    }

    /**
     * swap dimensions and device memory with another parallel exclusive prefix sum
     */
    void swap(prefix_sum& p)
    {
	std::swap(count, p.count);
	std::swap(threads, p.threads);
	blocks.swap(p.blocks);
	g_sum.swap(p.g_sum);
    }

    /**
     * in-place parallel exclusive prefix sum
     */
    void operator()(cuda::vector<T>& g_array, cuda::stream& stream)
    {
	assert(g_array.size() == count);

	// compute blockwise partial prefix sums
	cuda::configure(blocks[0], threads, gpu::prefix_sum::boff(2 * threads * sizeof(T)), stream);
	gpu::prefix_sum::grid_prefix_sum(g_array, g_array, g_sum[0], count);

	for (uint i = 1; i < blocks.size(); ++i) {
	    cuda::configure(blocks[i], threads, gpu::prefix_sum::boff(2 * threads * sizeof(T)), stream);
	    gpu::prefix_sum::grid_prefix_sum(g_sum[i - 1], g_sum[i - 1], g_sum[i], blocks[i - 1]);
	}

	// add block prefix sums to partial prefix sums
	if (blocks.size() > 1) {
	    for (uint i = blocks.size() - 2; i > 0; --i) {
		cuda::configure(g_sum[i].size(), threads, stream);
		gpu::prefix_sum::add_block_sums(g_sum[i - 1], g_sum[i - 1], g_sum[i], blocks[i - 1]);
	    }
	}
	if (blocks[0] > 1) {
	    cuda::configure(blocks[0], threads, stream);
	    gpu::prefix_sum::add_block_sums(g_array, g_array, g_sum[0], count);
	}
    }

private:
    uint count, threads;
    std::vector<cuda::vector<T> > g_sum;
    std::vector<uint> blocks;
};

} // namespace ljgpu

#endif /* ! LJGPU_ALGORITHM_PREFIX_SUM_HPP */
