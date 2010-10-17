/*
 * Copyright © 2008-2010  Peter Colberg
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

#include <halmd/mdsim/gpu/particle_kernel.cuh>
#include <halmd/mdsim/gpu/sampler/trajectory_kernel.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/utility/gpu/thread.cuh>
#include <halmd/utility/gpu/variant.cuh>

using namespace halmd::mdsim::gpu::particle_kernel;
using namespace halmd::utility::gpu;

namespace halmd
{
namespace mdsim { namespace gpu { namespace sampler
{
namespace trajectory_kernel
{

/** positions, types */
texture<float4> r_;
/** minimum image vectors */
texture<variant<map<pair<int_<3>, float4>, pair<int_<2>, float2> > > > image_;
/** velocities, tags */
texture<float4> v_;
/** cubic box edgle length */
__constant__ variant<map<pair<int_<3>, float3>, pair<int_<2>, float2> > > box_length_;

/**
 * sample trajectory for all particle of a single species
 */
template <typename vector_type, typename T>
__global__ void sample(unsigned int const* g_index, T* g_or, T* g_ov)
{
    enum { dimension = vector_type::static_size };
    // permutation index
    uint const j = g_index[GTID];
    // fetch particle from texture caches
    unsigned int tag, type;
    vector_type r, v;
    tie(r, type) = untagged<vector_type>(tex1Dfetch(r_, j));
    tie(v, tag) = untagged<vector_type>(tex1Dfetch(v_, j));
    vector_type image = tex1Dfetch(get<dimension>(image_), j);
    vector_type L = get<dimension>(box_length_);
    // store particle in global memory
    g_or[GTID] = r + element_prod(L, image);
    g_ov[GTID] = v;
}

} // namespace trajectory_kernel

template <int dimension>
trajectory_wrapper<dimension> const trajectory_wrapper<dimension>::kernel = {
    trajectory_kernel::r_
  , get<dimension>(trajectory_kernel::image_)
  , trajectory_kernel::v_
  , get<dimension>(trajectory_kernel::box_length_)
  , trajectory_kernel::sample<fixed_vector<float, dimension> >
};

template class trajectory_wrapper<3>;
template class trajectory_wrapper<2>;

}}} // namespace mdsim::gpu::sampler

} // namespace halmd
