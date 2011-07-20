/*
 * Copyright © 2011  Felix Höfling
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

#include <halmd/mdsim/gpu/box_kernel.cuh>
#include <test/unit/mdsim/box_kernel.hpp>

// thread ID within block
#define TID     threadIdx.x
// number of threads per block
#define TDIM    blockDim.x
// block ID within grid
#define BID     (blockIdx.y * gridDim.x + blockIdx.x)
// number of blocks within grid
#define BDIM    (gridDim.y * gridDim.x)
// thread ID within grid
#define GTID    (BID * TDIM + TID)
// number of threads per grid
#define GTDIM   (BDIM * TDIM)

namespace box {

using namespace halmd::mdsim::gpu;

/**
 * fill float array with a given value using a loop
 */
template <typename coalesced_vector_type, typename vector_type>
__global__ void reduce_periodic(coalesced_vector_type* g_r, coalesced_vector_type* g_reduced, const vector_type length)
{
    const unsigned int i = GTID;

    vector_type r = g_r[i];
    vector_type image = box_kernel::reduce_periodic(r, length);
    g_reduced[i] = r;

    box_kernel::extend_periodic(r, image, length);
    g_r[i] = r;
}

} // namespace box

template <int dimension, typename float_type>
box_kernel_wrapper<dimension, float_type> box_kernel_wrapper<dimension, float_type>::kernel = {
    box::reduce_periodic
};

template class box_kernel_wrapper<2, float>;
template class box_kernel_wrapper<3, float>;
