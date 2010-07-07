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

#ifndef HALMD_ALGORITHM_GPU_SCAN_HPP
#define HALMD_ALGORITHM_GPU_SCAN_HPP

#include <stdexcept>
#include <vector>

#include <halmd/algorithm/gpu/scan_kernel.hpp>

namespace halmd
{
namespace algorithm { namespace gpu
{

/*
 * Parallel exclusive prefix sum
 */
template <typename T>
class scan
{
public:
    /**
     * allocate parallel exclusive prefix sum for given element count
     */
    scan(uint const& count, uint const& threads)
      : count(count)
      , threads(threads)
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
     * in-place parallel exclusive prefix sum
     */
    void operator()(cuda::vector<T>& g_array)
    {
        assert(g_array.size() == count);

        // compute blockwise partial prefix sums
        cuda::configure(blocks[0], threads, boff(2 * threads * sizeof(T)));
        get_scan_kernel<T>().grid_prefix_sum(g_array, g_array, g_sum[0], count);

        for (uint i = 1; i < blocks.size(); ++i) {
            cuda::configure(blocks[i], threads, boff(2 * threads * sizeof(T)));
            get_scan_kernel<T>().grid_prefix_sum(g_sum[i - 1], g_sum[i - 1], g_sum[i], blocks[i - 1]);
        }

        // add block prefix sums to partial prefix sums
        if (blocks.size() > 1) {
            for (uint i = blocks.size() - 2; i > 0; --i) {
                cuda::configure(g_sum[i].size(), threads);
                get_scan_kernel<T>().add_block_sums(g_sum[i - 1], g_sum[i - 1], g_sum[i], blocks[i - 1]);
            }
        }
        if (blocks[0] > 1) {
            cuda::configure(blocks[0], threads);
            get_scan_kernel<T>().add_block_sums(g_array, g_array, g_sum[0], count);
        }
    }

private:
    uint count, threads;
    std::vector<cuda::vector<T> > g_sum;
    std::vector<uint> blocks;
};

}} // namespace algorithm::gpu

} // namespace halmd

#endif /* ! HALMD_ALGORITHM_GPU_SCAN_HPP */
