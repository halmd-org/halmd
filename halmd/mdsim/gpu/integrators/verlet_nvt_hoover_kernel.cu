/*
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
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

// this file is identical to verlet_kernel.cu except for the
// rescaling of velocities in integrate()

#include <boost/mpl/if.hpp>

#include <halmd/algorithm/gpu/reduce_kernel.cuh>
#include <halmd/mdsim/gpu/integrators/verlet_kernel.cuh>
#include <halmd/mdsim/gpu/integrators/verlet_nvt_hoover_kernel.hpp>
#include <halmd/mdsim/gpu/particle_kernel.cuh>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>
#include <halmd/utility/gpu/thread.cuh>

using namespace boost::mpl;
using namespace halmd::algorithm::gpu;
using namespace halmd::mdsim::gpu::particle_kernel;

namespace halmd
{
namespace mdsim { namespace gpu { namespace integrators
{
namespace verlet_nvt_hoover_kernel
{

/** integration time-step */
static __constant__ float timestep_;

/**
 * First leapfrog half-step of velocity-Verlet algorithm
 */
template <
    typename float_type
  , typename vector_type
  , typename vector_type_
  , typename gpu_vector_type
>
__global__ void _integrate(
    float4* g_r
  , gpu_vector_type* g_image
  , float4* g_v
  , gpu_vector_type const* g_f
  , float_type scale
  , vector_type_ box_length
)
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

    v *= scale;          //< rescale velocity according to Nosé-Hoover dynamics
    verlet_kernel::integrate(r, image, v, f, timestep_, box_length);

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
__global__ void _finalize(float4* g_v, gpu_vector_type const* g_f)
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

    verlet_kernel::finalize(v, f, timestep_);

#ifdef USE_VERLET_DSFUN
    tie(g_v[i], g_v[i + threads]) = tagged(v, tag);
#else
    g_v[i] = tagged(v, tag);
#endif
}

/**
 * rescale velocities
 */
template <typename float_type, typename vector_type>
__global__ void rescale(float4* g_v, float_type scale)
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

    v *= scale;

#ifdef USE_VERLET_DSFUN
    tie(g_v[i], g_v[i + threads]) = tagged(v, tag);
#else
    g_v[i] = tagged(v, tag);
#endif
}

} // namespace verlet_nvt_hoover_kernel

template <int dimension, typename float_type>
verlet_nvt_hoover_wrapper<dimension, float_type> const
verlet_nvt_hoover_wrapper<dimension, float_type>::kernel = {
    verlet_nvt_hoover_kernel::timestep_
  , verlet_nvt_hoover_kernel::_integrate<
        float_type
      , fixed_vector<float_type, dimension>
      , fixed_vector<float, dimension>
    >
  , verlet_nvt_hoover_kernel::_finalize<
        fixed_vector<float_type, dimension>
      , fixed_vector<float, dimension>
    >
  , verlet_nvt_hoover_kernel::rescale<float_type, fixed_vector<float_type, dimension> >
};

#ifdef USE_VERLET_DSFUN
template class verlet_nvt_hoover_wrapper<3, dsfloat>;
template class verlet_nvt_hoover_wrapper<2, dsfloat>;
#else
template class verlet_nvt_hoover_wrapper<3, float>;
template class verlet_nvt_hoover_wrapper<2, float>;
#endif

}}} // namespace mdsim::gpu::integrators

template class reduce_wrapper<
    sum_                        // reduce_transform
  , fixed_vector<float, 3>      // input_type
  , float4                      // coalesced_input_type
  , dsfloat                     // output_type
  , dsfloat                     // coalesced_output_type
  , square_                     // input_transform
>;

template class reduce_wrapper<
    sum_                        // reduce_transform
  , fixed_vector<float, 2>      // input_type
  , float4                      // coalesced_input_type
  , dsfloat                     // output_type
  , dsfloat                     // coalesced_output_type
  , square_                     // input_transform
>;

} // namespace halmd
