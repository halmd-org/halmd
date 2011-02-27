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
#include <halmd/mdsim/gpu/positions/trajectory_kernel.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/utility/gpu/thread.cuh>
#include <halmd/utility/gpu/variant.cuh>

using namespace halmd::mdsim::gpu::particle_kernel;
using namespace halmd::utility::gpu;

/** cuboid box edge length */
static __constant__ variant<map<pair<int_<3>, float3>, pair<int_<2>, float2> > > box_length_;

namespace halmd
{
namespace mdsim { namespace gpu { namespace positions
{
namespace trajectory_kernel
{

/**
 * shift particle positions to range (-L/2, L/2)
 */
template <typename vector_type>
__global__ void reduce_periodic(float4* g_r)
{
    enum { dimension = vector_type::static_size };

    vector_type r;
    unsigned int type;
    tie(r, type) = untagged<vector_type>(g_r[GTID]);

    vector_type box_length = get<dimension>(box_length_);

    box_kernel::reduce_periodic(r, box_length);

    g_r[GTID] = tagged(r, type);
}

} // namespace trajectory_kernel

template <int dimension>
trajectory_wrapper<dimension> const trajectory_wrapper<dimension>::kernel = {
    get<dimension>(box_length_)
  , trajectory_kernel::reduce_periodic<fixed_vector<float, dimension> >
};

template class trajectory_wrapper<3>;
template class trajectory_wrapper<2>;

}}} // namespace mdsim::gpu::positions

} // namespace halmd
