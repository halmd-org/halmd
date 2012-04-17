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
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>
#include <halmd/utility/gpu/thread.cuh>

using namespace boost::mpl;
using namespace halmd::mdsim::gpu::particle_kernel;

namespace halmd {
namespace mdsim {
namespace gpu {
namespace integrators {
namespace verlet_kernel {

/** integration time-step */
static __constant__ float timestep_;

/**
 * First leapfrog half-step of velocity-Verlet algorithm
 */
template <
    typename vector_type
  , typename vector_type_
  , typename gpu_vector_type
>
__global__ void _integrate(
    float4* g_r
  , gpu_vector_type* g_image
  , float4* g_v
  , gpu_vector_type const* g_f
  , float const* g_mass
  , unsigned int ntype
  , vector_type_ box_length
)
{
    extern __shared__ float s_mass[];
    if (TID < ntype) {
        s_mass[TID] = g_mass[TID];
    }
    __syncthreads();

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
    float mass = s_mass[type];

    integrate(r, image, v, f, mass, timestep_, box_length);

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
    float4 const* g_r
  , float4* g_v
  , gpu_vector_type const* g_f
  , float const* g_mass
  , unsigned int ntype
)
{
    extern __shared__ float s_mass[];
    if (TID < ntype) {
        s_mass[TID] = g_mass[TID];
    }
    __syncthreads();

    unsigned int const i = GTID;
    unsigned int const threads = GTDIM;
    unsigned int tag, type;
    vector_type v;
    vector_type_ _;
    tie(_, type) = untagged<vector_type_>(g_r[i]);
#ifdef USE_VERLET_DSFUN
    tie(v, tag) = untagged<vector_type>(g_v[i], g_v[i + threads]);
#else
    tie(v, tag) = untagged<vector_type>(g_v[i]);
#endif
    vector_type_ f = g_f[i];
    float mass = s_mass[type];

    finalize(v, f, mass, timestep_);

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
#ifdef USE_VERLET_DSFUN
  , verlet_kernel::_integrate<fixed_vector<dsfloat, dimension>, fixed_vector<float, dimension> >
  , verlet_kernel::_finalize<fixed_vector<dsfloat, dimension>, fixed_vector<float, dimension> >
#else
  , verlet_kernel::_integrate<fixed_vector<float, dimension>, fixed_vector<float, dimension> >
  , verlet_kernel::_finalize<fixed_vector<float, dimension>, fixed_vector<float, dimension> >
#endif
};

template class verlet_wrapper<3>;
template class verlet_wrapper<2>;

} // namespace mdsim
} // namespace gpu
} // namespace integrators
} // namespace halmd
