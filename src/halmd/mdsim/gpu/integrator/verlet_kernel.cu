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
#include <halmd/mdsim/gpu/integrator/verlet_kernel.cuh>
#include <halmd/mdsim/gpu/particle_kernel.cuh>
#include <halmd/numeric/gpu/blas/vector.cuh>

using namespace boost::mpl;
using namespace halmd::mdsim::gpu::particle_kernel;
using namespace halmd::numeric::gpu::blas;

namespace halmd { namespace mdsim { namespace gpu { namespace integrator
{

namespace verlet_kernel
{

/** integration time-step */
__constant__ float timestep_;

template <size_t N>
struct dim_
{
    /** cubic box edgle length */
    static __constant__ typename if_c<N == 3, float3, float2>::type box_length;
};

// explicit instantiation
template class dim_<3>;
template class dim_<2>;

/**
 * First leapfrog half-step of velocity-Verlet algorithm
 */
template <typename vector_type, typename vector_type_, typename gpu_vector_type>
__global__ void _integrate(
  float4* g_r,
  gpu_vector_type* g_image,
  float4* g_v,
  gpu_vector_type const* g_f)
{
    unsigned int const i = GTID;
    unsigned int const threads = GTDIM;
    unsigned int type, tag;
#ifdef USE_VERLET_DSFUN
    vector_type r      ( untagged<vector_type_>(g_r[i], type),
                         untagged<vector_type_>(g_r[i + threads]));
    vector_type v      ( untagged<vector_type_>(g_v[i], tag),
                         untagged<vector_type_>(g_v[i + threads]));
#else
    vector_type_ r     = untagged<vector_type_>(g_r[i], type);
    vector_type_ v     = untagged<vector_type_>(g_v[i], tag);
#endif
    vector_type_ image = g_image[i];
    vector_type_ f     = g_f[i];
    vector_type_ L     = dim_<vector_type::static_size>::box_length;

    integrate(r, image, v, f, timestep_, L);

#ifdef USE_VERLET_DSFUN
    g_r[i]             = tagged(dsfloat_hi(r), type);
    g_r[i + threads]   = tagged(dsfloat_lo(r), 0);
    g_v[i]             = tagged(dsfloat_hi(v), tag);
    g_v[i + threads]   = tagged(dsfloat_lo(v), 0);
#else
    g_r[i]             = tagged(r, type);
    g_v[i]             = tagged(v, tag);
#endif
    g_image[i]         = image;
}

/**
 * Second leapfrog half-step of velocity-Verlet algorithm
 */
template <typename vector_type, typename vector_type_, typename gpu_vector_type>
__global__ void _finalize(
  float4* g_v,
  gpu_vector_type const* g_f)
{
    unsigned int const i = GTID;
    unsigned int const threads = GTDIM;
    unsigned int tag;
#ifdef USE_VERLET_DSFUN
    vector_type v      ( untagged<vector_type_>(g_v[i], tag),
                         untagged<vector_type_>(g_v[i + threads]));
#else
    vector_type_ v     = untagged<vector_type_>(g_v[i], tag);
#endif
    vector_type_ f     = g_f[i];

    finalize(v, f, timestep_);

#ifdef USE_VERLET_DSFUN
    g_v[i]             = tagged(dsfloat_hi(v), tag);
    g_v[i + threads]   = tagged(dsfloat_lo(v), 0);
#else
    g_v[i]             = tagged(v, tag);
#endif
}

} // namespace verlet_kernel

}}}} //namespace halmd::mdsim::gpu::integrator
