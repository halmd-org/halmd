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

#ifndef LJGPU_MDSIM_VELOCITY_HPP
#define LJGPU_MDSIM_VELOCITY_HPP

#include <algorithm>
#include <boost/array.hpp>
#include <cuda_wrapper.hpp>
#include <ljgpu/mdsim/gpu/velocity.hpp>

namespace ljgpu
{

/*
 * Parallel reduction
 */
template <typename gpu_output_type, typename output_type>
class squared_velocity_sum
{
public:
    enum { BLOCKS = gpu::velocity::BLOCKS };
    enum { THREADS = gpu::velocity::THREADS };

    /**
     * stream parallel reduction kernel
     */
    template <typename T>
    void operator()(T const& g_in, unsigned int offset, cuda::stream& stream)
    {
	g_block_sum.resize(BLOCKS);
	h_block_sum.resize(BLOCKS);
	cuda::configure(BLOCKS, THREADS, stream);
	gpu::velocity::sum(g_in, g_block_sum, g_in.size(), offset);
	cuda::copy(g_block_sum, h_block_sum, stream);
    }

    /**
     * returns sum after CUDA stream synchronisation
     */
    output_type value() const
    {
	return std::accumulate(h_block_sum.begin(), h_block_sum.end(), output_type(0));
    }

private:
    cuda::vector<gpu_output_type> g_block_sum;
    cuda::host::vector<gpu_output_type> h_block_sum;
};

} // namespace ljgpu

#endif /* ! LJGPU_MDSIM_VELOCITY_HPP */
