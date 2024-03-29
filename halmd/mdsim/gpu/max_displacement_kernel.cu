/*
 * Copyright © 2008-2011  Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <halmd/algorithm/gpu/reduction.cuh>
#include <halmd/mdsim/gpu/box_kernel.cuh>
#include <halmd/mdsim/gpu/max_displacement_kernel.hpp>
#include <halmd/utility/gpu/thread.cuh>

using namespace halmd::algorithm::gpu;

namespace halmd {
namespace mdsim {
namespace gpu {
namespace max_displacement_kernel {

/**
 * maximum squared particle displacement
 */
template <typename vector_type>
__global__ void displacement(
    float4 const* g_r
  , float4 const* g_r0
  , typename vector_type::value_type* g_rr
  , unsigned int npart
  , vector_type box_length
)
{
    typedef typename vector_type::value_type float_type;
    enum { dimension = vector_type::static_size };

    float_type rr = 0;

    for (uint i = GTID; i < npart; i += GTDIM) {
        vector_type r;
        unsigned int type;
        tie(r, type) <<= g_r[i];
        vector_type r0;
        tie(r0, type) <<= g_r0[i];
        r -= r0;
        box_kernel::reduce_periodic(r, box_length);
        rr = max(rr, inner_prod(r, r));
    }

    // reduce values for all threads in block with the maximum function
    reduce<max_>(rr);

    if (TID < 1) {
        // store block reduced value in global memory
        g_rr[blockIdx.x] = rr;
    }
}

} // namespace max_displacement_kernel

template <int dimension>
max_displacement_wrapper<dimension> max_displacement_wrapper<dimension>::kernel = {
    max_displacement_kernel::displacement<fixed_vector<float, dimension>>
};

template class max_displacement_wrapper<3>;
template class max_displacement_wrapper<2>;

} //namespace gpu
} //namespace mdsim
} //namespace halmd
