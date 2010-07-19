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

#include <halmd/mdsim/gpu/integrators/verlet_kernel.cuh>
#include <halmd/mdsim/gpu/integrators/verlet_kernel.hpp>
#include <halmd/mdsim/gpu/particle_kernel.cuh>
#include <halmd/numeric/gpu/blas/dsfloat.cuh>
#include <halmd/numeric/gpu/blas/vector.cuh>
#include <halmd/utility/gpu/thread.cuh>
#include <halmd/utility/gpu/variant.cuh>

using namespace boost::mpl;
using namespace halmd::mdsim::gpu::particle_kernel;
using namespace halmd::numeric::gpu::blas;
using namespace halmd::utility::gpu;

namespace halmd
{
namespace mdsim { namespace gpu { namespace integrators
{
namespace verlet_kernel
{

/** integration time-step */
static __constant__ float timestep_;
/** cuboid box edge length */
static __constant__ variant<map<pair<int_<3>, float3>, pair<int_<2>, float2> > > box_length_;

/**
 * First leapfrog half-step of velocity-Verlet algorithm
 */
template <
    typename vector_type
  , typename vector_type_
  , typename gpu_vector_type
>
__global__ void _integrate(
  float4* g_r,
  gpu_vector_type* g_image,
  float4* g_v,
  gpu_vector_type const* g_f)
{
    unsigned int const i = GTID;
    unsigned int const threads = GTDIM;
    unsigned int type, tag;
    vector_type r, v;
#ifdef USE_VERLET_DSFUN
    tie(r, type) = untagged<vector_type>(g_r[i], g_r[i + threads]);
    tie(v, tag) = untagged<vector_type>(g_v[i], g_v[i + threads]);
#else
    tie(r, type) = untagged<vector_type>(g_r[i]);
    tie(v, tag) = untagged<vector_type>(g_v[i]);
#endif
    vector_type_ image = g_image[i];
    vector_type_ f = g_f[i];
    vector_type_ L = get<vector_type::static_size>(box_length_);

    integrate(r, image, v, f, timestep_, L);

#ifdef USE_VERLET_DSFUN
    tie(g_r[i], g_r[i + threads]) = tagged(r, type);
    tie(g_v[i], g_v[i + threads]) = tagged(v, tag);
#else
    g_r[i] = tagged(r, type);
    g_v[i] = tagged(v, tag);
#endif
    g_image[i] = image;
}

/**
 * Second leapfrog half-step of velocity-Verlet algorithm
 */
template <
    typename vector_type
  , typename vector_type_
  , typename gpu_vector_type
>
__global__ void _finalize(
  float4* g_v,
  gpu_vector_type const* g_f)
{
    unsigned int const i = GTID;
    unsigned int const threads = GTDIM;
    unsigned int tag;
    vector_type v;
#ifdef USE_VERLET_DSFUN
    tie(v, tag) = untagged<vector_type>(g_v[i], g_v[i + threads]);
#else
    tie(v, tag) = untagged<vector_type>(g_v[i]);
#endif
    vector_type_ f = g_f[i];

    finalize(v, f, timestep_);

#ifdef USE_VERLET_DSFUN
    tie(g_v[i], g_v[i + threads]) = tagged(v, tag);
#else
    g_v[i] = tagged(v, tag);
#endif
}

} // namespace verlet_kernel

template <int dimension>
verlet_wrapper<dimension> const verlet_wrapper<dimension>::wrapper = {
    verlet_kernel::timestep_
  , get<dimension>(verlet_kernel::box_length_)
#ifdef USE_VERLET_DSFUN
  , verlet_kernel::_integrate<vector<dsfloat, dimension>, vector<float, dimension> >
  , verlet_kernel::_finalize<vector<dsfloat, dimension>, vector<float, dimension> >
#else
  , verlet_kernel::_integrate<vector<float, dimension>, vector<float, dimension> >
  , verlet_kernel::_finalize<vector<float, dimension>, vector<float, dimension> >
#endif
};

template class verlet_wrapper<3>;
template class verlet_wrapper<2>;

}}} // namespace mdsim::gpu::integrators

} // namespace halmd
