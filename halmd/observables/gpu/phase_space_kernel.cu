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

#include <halmd/mdsim/gpu/box_kernel.cuh>
#include <halmd/mdsim/gpu/particle_kernel.cuh>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/observables/gpu/phase_space_kernel.hpp>
#include <halmd/utility/gpu/thread.cuh>
#include <halmd/utility/gpu/variant.cuh>

using namespace halmd::utility::gpu; //< variant, map, pair
using namespace halmd::mdsim::gpu; //< namespace box_kernel

namespace halmd {
namespace observables {
namespace gpu {
namespace phase_space_kernel {

/** positions, types */
texture<float4> r_;
/** minimum image vectors */
texture<variant<map<pair<int_<3>, float4>, pair<int_<2>, float2> > > > image_;
/** velocities, tags */
texture<float4> v_;
/** cubic box edgle length */
__constant__ variant<map<pair<int_<3>, float3>, pair<int_<2>, float2> > > box_length_;

/**
 * sample phase space for all particle of a single species
 */
template <typename vector_type, typename T>
__global__ void sample(unsigned int const* g_index, T* g_or, T* g_ov, unsigned int ntype)
{
    using mdsim::gpu::particle_kernel::untagged;

    enum { dimension = vector_type::static_size };

    if (GTID < ntype) {
        // permutation index
        uint const j = g_index[GTID];
        // fetch particle from texture caches
        unsigned int tag, type;
        vector_type r, v;
        tie(r, type) = untagged<vector_type>(tex1Dfetch(r_, j));
        tie(v, tag) = untagged<vector_type>(tex1Dfetch(v_, j));
        // extend particle positions in periodic box
        vector_type image = tex1Dfetch(get<dimension>(image_), j);
        vector_type L = get<dimension>(box_length_);
        box_kernel::extend_periodic(r, image, L);
        // store particle in global memory
        g_or[GTID] = r;
        g_ov[GTID] = v;
    }
}

} // namespace phase_space_kernel

template <int dimension>
phase_space_wrapper<dimension> const phase_space_wrapper<dimension>::kernel = {
    phase_space_kernel::r_
  , get<dimension>(phase_space_kernel::image_)
  , phase_space_kernel::v_
  , get<dimension>(phase_space_kernel::box_length_)
  , phase_space_kernel::sample<fixed_vector<float, dimension> >
};

template class phase_space_wrapper<3>;
template class phase_space_wrapper<2>;

} // namespace observables
} // namespace gpu
} // namespace halmd
