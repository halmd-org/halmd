/* Parallel radix sort
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

#ifndef MDSIM_GPU_RADIX_GLUE_HPP
#define MDSIM_GPU_RADIX_GLUE_HPP

#include <cuda_wrapper.hpp>

namespace mdsim { namespace gpu { namespace radix
{

enum {
    // number of threads in half-warp
    HALF_WARP_SIZE = 16,
    // radix bits per counting sort iteration
    RADIX = 8,
    // bucket count per partition
    BUCKET_SIZE = (1 << RADIX),
    // bucket count per thread
    BUCKETS_PER_THREAD = BUCKET_SIZE / HALF_WARP_SIZE,
};

extern cuda::function<
    void (uint const*, uint*, const uint, const uint)
    > histogram_keys;

extern cuda::function<
    void (uint const*, uint*, uint const*, uint*, uint const*, const uint, const uint),
    void (uint const*, uint*, int const*, int*, uint const*, const uint, const uint),
    void (uint const*, uint*, float2 const*, float2*, uint const*, const uint, const uint),
    void (uint const*, uint*, float4 const*, float4*, uint const*, const uint, const uint)
    > permute;

}}} // namespace mdsim::gpu::radix

#endif /* ! MDSIM_GPU_RADIX_GLUE_HPP */
