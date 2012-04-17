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

#include <halmd/mdsim/gpu/box_kernel.cuh>
#include <halmd/mdsim/gpu/particle_kernel.cuh>
#include <halmd/mdsim/gpu/positions/phase_space_kernel.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/utility/gpu/thread.cuh>

using namespace halmd::mdsim::gpu::particle_kernel;

namespace halmd {
namespace mdsim {
namespace gpu {
namespace positions {
namespace phase_space_kernel {

/**
 * shift particle positions to range (-L/2, L/2)
 *
 * FIXME move to box_kernel.cuh
 */
template <typename vector_type, typename coalesced_vector_type>
__global__ void reduce_periodic(
    float4* g_r
  , coalesced_vector_type* g_image
  , vector_type box_length
)
{
    enum { dimension = vector_type::static_size };

    vector_type r;
    unsigned int type;
    tie(r, type) = untagged<vector_type>(g_r[GTID]);

    vector_type image = box_kernel::reduce_periodic(r, box_length);

    g_image[GTID] = image;
    g_r[GTID] = tagged(r, type);
}

} // namespace phase_space_kernel

template <int dimension>
phase_space_wrapper<dimension> const phase_space_wrapper<dimension>::kernel = {
    phase_space_kernel::reduce_periodic<fixed_vector<float, dimension> >
};

template class phase_space_wrapper<3>;
template class phase_space_wrapper<2>;

} // namespace mdsim
} // namespace gpu
} // namespace positions
} // namespace halmd
