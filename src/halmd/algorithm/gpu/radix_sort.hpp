/* Parallel radix sort
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

#ifndef HALMD_ALGORITHM_GPU_RADIX_SORT_HPP
#define HALMD_ALGORITHM_GPU_RADIX_SORT_HPP

#include <cuda_wrapper.hpp>

namespace halmd { namespace gpu { namespace radix_sort
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

extern cuda::function<void (uint const*, uint*, uint, uint)> histogram_keys;
extern cuda::function<void (uint const*, uint*, int const*, int*, uint const*, uint, uint),
                      void (uint const*, uint*, uint const*, uint*, uint const*, uint, uint),
                      void (uint const*, uint*, float const*, float*, uint const*, uint, uint),
                      void (uint const*, uint*, float2 const*, float2*, uint const*, uint, uint),
                      void (uint const*, uint*, float4 const*, float4*, uint const*, uint, uint)> permute;

}}} // namespace halmd::gpu::radix_sort

#endif /* ! HALMD_ALGORITHM_GPU_RADIX_SORT_HPP */
