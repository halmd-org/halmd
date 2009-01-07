/* Parallel exclusive prefix sum
 *
 * Copyright (C) 2008  Peter Colberg
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

#ifndef MDSIM_GPU_SCAN_GLUE_HPP
#define MDSIM_GPU_SCAN_GLUE_HPP

#include <cuda_wrapper.hpp>

namespace mdsim { namespace gpu { namespace scan
{

enum {
    // number of shared memory banks
    SHMEM_BANKS = 16,
};

/**
 * returns array index with offset for bank conflict free shared memory access
 */
__device__ __host__ inline uint boff(uint const& i)
{
    return i + i / SHMEM_BANKS;
}

template <typename T>
struct __scan
{
    static cuda::function<void (T const*, T*, T*, const uint)> block_prefix_sum;
    static cuda::function<void (T const*, T*, const uint)> prefix_sum;
    static cuda::function<void (T const*, T*, T const*, const uint)> add_block_sums;
};

template <typename T0, typename T1, typename T2, typename T3>
inline void block_prefix_sum(T0 const& t0, T1& t1, T2& t2, T3 const& t3)
{
    __scan<typename T0::value_type>::block_prefix_sum(t0, t1, t2, t3);
}

template <typename T0, typename T1, typename T2>
inline void prefix_sum(T0 const& t0, T1& t1, T2 const& t2)
{
    __scan<typename T0::value_type>::prefix_sum(t0, t1, t2);
}

template <typename T0, typename T1, typename T2, typename T3>
inline void add_block_sums(T0 const& t0, T1& t1, T2 const& t2, T3 const& t3)
{
    __scan<typename T0::value_type>::add_block_sums(t0, t1, t2, t3);
}

}}} // namespace mdsim::gpu::scan

#endif /* ! MDSIM_GPU_SCAN_GLUE_HPP */
