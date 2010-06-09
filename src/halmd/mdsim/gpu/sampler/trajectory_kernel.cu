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

#include <boost/mpl/if.hpp>

#include <halmd/algorithm/gpu/base.cuh>
#include <halmd/mdsim/gpu/particle_kernel.cuh>
#include <halmd/numeric/gpu/blas/vector.cuh>

using namespace boost::mpl;
using namespace halmd::mdsim::gpu::particle_kernel;
using namespace halmd::numeric::gpu::blas;

namespace halmd { namespace mdsim { namespace gpu { namespace sampler { namespace trajectory_kernel
{

template <size_t N>
struct dim_
{
    typedef typename if_c<N == 3, float4, float2>::type coalesced_vector_type;
    typedef typename if_c<N == 3, float3, float2>::type vector_type;

    /** positions, types */
    static texture<float4, 1, cudaReadModeElementType> r;
    /** minimum image vectors */
    static texture<coalesced_vector_type, 1, cudaReadModeElementType> image;
    /** velocities, tags */
    static texture<float4, 1, cudaReadModeElementType> v;
    /** cubic box edgle length */
    static __constant__ vector_type box_length;
};

// explicit instantiation
template class dim_<3>;
template class dim_<2>;

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
    vector_type r = untagged<vector_type >(tex1Dfetch(dim_<dimension>::r, j));
    vector_type image = tex1Dfetch(dim_<dimension>::image, j);
    vector_type L = dim_<dimension>::box_length;
    vector_type v = untagged<vector_type >(tex1Dfetch(dim_<dimension>::v, j));
    // store particle in global memory
    g_or[GTID] = r + element_prod(L, image);
    g_ov[GTID] = v;
}

}}}}} // namespace halmd::mdsim::gpu::sampler::trajectory_kernel
