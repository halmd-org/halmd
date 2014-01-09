/*
 * Copyright © 2014 Felix Höfling
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

#include <halmd/observables/modulation.hpp>
#include <test/unit/observables/modulation_kernel.hpp>

using namespace halmd;

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

namespace modulation_kernel {

template <typename coalesced_vector_type, typename modulation_type>
__global__ void compute(
    coalesced_vector_type const* g_r
  , float* g_out
  , unsigned int const size
  , modulation_type const modulation
)
{
    typedef typename modulation_type::argument_type vector_type;
    unsigned int const i = GTID;

    if (i < size) {
        vector_type r = g_r[i];
        g_out[i] = modulation(r);
    }
}

} // namespace modulation_kernel

template <int dimension, typename modulation_type>
modulation_wrapper<dimension, modulation_type> const
modulation_wrapper<dimension, modulation_type>::kernel = {
    modulation_kernel::compute
};

using namespace halmd::observables::modulation;

template class modulation_wrapper<3, exponential<3, float> >;
template class modulation_wrapper<2, exponential<2, float> >;
template class modulation_wrapper<3, exponential_slab<3, float> >;
template class modulation_wrapper<2, exponential_slab<2, float> >;
template class modulation_wrapper<3, catenary<3, float> >;
template class modulation_wrapper<2, catenary<2, float> >;
