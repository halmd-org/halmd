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

#ifndef HALMD_MDSIM_GPU_INTEGRATOR_VERLET_NVT_HOOVER_KERNEL_HPP
#define HALMD_MDSIM_GPU_INTEGRATOR_VERLET_NVT_HOOVER_KERNEL_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <halmd/mdsim/type_traits.hpp>

namespace halmd
{
namespace mdsim { namespace gpu { namespace integrators
{

template <int dimension, typename float_type>
struct verlet_nvt_hoover_wrapper
{
    typedef typename type_traits<dimension, float>::gpu::coalesced_vector_type coalesced_vector_type;
    typedef typename type_traits<dimension, float>::gpu::vector_type vector_type;

    /** integration time-step */
    cuda::symbol<float> timestep;
    /** cubic box edge length */
    cuda::symbol<vector_type> box_length;
    /** first leapfrog half-step of velocity-Verlet algorithm */
    cuda::function <void (float4*, coalesced_vector_type*, float4*, coalesced_vector_type const*, float_type)> integrate;
    /** second leapfrog half-step of velocity-Verlet algorithm */
    cuda::function <void (float4*, coalesced_vector_type const*)> finalize;
    /** rescaling of velocities */
    cuda::function <void (float4*, float_type)> rescale;

    static verlet_nvt_hoover_wrapper const kernel;
};

}}} // namespace mdsim::gpu::integrators

} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_INTEGRATOR_VERLET_NVT_HOOVER_KERNEL_HPP */
