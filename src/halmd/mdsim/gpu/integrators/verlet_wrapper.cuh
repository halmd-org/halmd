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

#ifndef HALMD_MDSIM_GPU_INTEGRATOR_VERLET_WRAPPER_CUH
#define HALMD_MDSIM_GPU_INTEGRATOR_VERLET_WRAPPER_CUH

#include <boost/mpl/if.hpp>

#include <cuda_wrapper.hpp>

namespace halmd { namespace mdsim { namespace gpu { namespace integrators
{

template <size_t N>
struct verlet_wrapper
{
    typedef typename boost::mpl::if_c<N == 3, float4, float2>::type coalesced_vector_type;
    typedef typename boost::mpl::if_c<N == 3, float3, float2>::type vector_type;

    /** integration time-step */
    static cuda::symbol<float> timestep;
    /** cubic box edgle length */
    static cuda::symbol<vector_type> length;
    /** first leapfrog half-step of velocity-Verlet algorithm */
    static cuda::function <void (float4*, coalesced_vector_type*, float4*, coalesced_vector_type const*)> integrate;
    /** second leapfrog half-step of velocity-Verlet algorithm */
    static cuda::function <void (float4*, coalesced_vector_type const*)> finalize;
};

}}}} // namespace halmd::mdsim::gpu::integrators

#endif /* ! HALMD_MDSIM_GPU_INTEGRATOR_VERLET_WRAPPER_CUH */
