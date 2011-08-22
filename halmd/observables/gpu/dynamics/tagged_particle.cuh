/*
 * Copyright © 2008-2011  Peter Colberg
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

#ifndef HALMD_OBSERVABLES_GPU_DYNAMICS_TAGGED_PARTICLE_CUH
#define HALMD_OBSERVABLES_GPU_DYNAMICS_TAGGED_PARTICLE_CUH

#include <halmd/algorithm/gpu/reduce_kernel.cuh>

namespace halmd {
namespace observables {
namespace gpu {
namespace dynamics {

template <
    typename correlation_function
  , unsigned int threads
  , typename coalesced_vector_type
  , typename accumulator_type
>
__global__ void accumulate(
    coalesced_vector_type const* g_r1
  , coalesced_vector_type const* g_r2
  , unsigned int size
  , accumulator_type* g_acc
)
{
    __shared__ accumulator_type s_acc[threads];
    accumulator_type acc;

    // load values from global device memory
    correlation_function correlate;
    for (unsigned int i = GTID; i < size; i += GTDIM) {
        acc(correlate(g_r1[i], g_r2[i]));
    }
    // reduced value for this thread
    s_acc[TID] = acc;
    __syncthreads();

    // compute reduced value for all threads in block
    algorithm::gpu::reduce<threads / 2, algorithm::gpu::accumulate_>(acc, acc, s_acc, s_acc);

    if (TID < 1) {
        // store block reduced value in global memory
        g_acc[blockIdx.x] = acc;
    }
}

} // namespace dynamics
} // namespace gpu
} // namespace observables
} // namespace halmd

#endif /* ! HALMD_OBSERVABLES_GPU_DYNAMICS_TAGGED_PARTICLE_CUH */
