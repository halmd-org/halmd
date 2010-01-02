/* Parallel reduction kernel
 *
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

#ifndef HALMD_MDSIM_VIRIAL_HPP
#define HALMD_MDSIM_VIRIAL_HPP

#include <algorithm>
#include <boost/array.hpp>

#include <cuda_wrapper.hpp>
#include <halmd/mdsim/gpu/virial.hpp>

namespace halmd
{

/*
 * Parallel reduction
 */
template <typename gpu_output_type, typename output_type>
class virial_sum
{
public:
    enum { BLOCKS = gpu::virial::BLOCKS };
    enum { THREADS = gpu::virial::THREADS };

    /**
     * stream parallel reduction kernel
     */
    template <typename T>
    void operator()(T const& g_in, T const& g_v, cuda::stream& stream)
    {
        g_block_sum.resize(BLOCKS);
        h_block_sum.resize(BLOCKS);
        cuda::configure(BLOCKS, THREADS, stream);
        gpu::virial::sum(g_in, g_v, g_block_sum, g_in.size());
        cuda::copy(g_block_sum, h_block_sum, stream);
    }

    template <typename T, typename U>
    void operator()(T const& g_in, T const& g_v, U const& g_tag,
                    boost::array<uint, 2> const& mpart, cuda::stream& stream)
    {
        g_block_sum.resize(2 * BLOCKS);
        h_block_sum.resize(2 * BLOCKS);
        cuda::configure(BLOCKS, THREADS, stream);
        gpu::virial::sum(g_in, g_v, g_tag, g_block_sum, g_in.size(), mpart.front());
        cuda::copy(g_block_sum, h_block_sum, stream);
    }

    /**
     * returns sum after CUDA stream synchronisation
     */
    std::vector<output_type> value() const
    {
        typedef typename cuda::host::vector<gpu_output_type>::const_iterator iterator;
        typedef typename output_type::iterator output_iterator;
        typedef typename output_type::value_type value_type;

        std::vector<output_type> v;
        for (iterator sum = h_block_sum.begin(); sum != h_block_sum.end(); sum += BLOCKS) {
            v.push_back(std::accumulate(sum, sum + BLOCKS, output_type(0)));
        }
        return v;
    }

private:
    cuda::vector<gpu_output_type> g_block_sum;
    cuda::host::vector<gpu_output_type> h_block_sum;
};

} // namespace halmd

#endif /* ! HALMD_MDSIM_VIRIAL_HPP */
