/* Lennard-Jones fluid kernel
 *
 * Copyright (C) 2008  Peter Colberg
 *
 * This program is free software: you can redistribute it and/or modify
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

#ifndef MDSIM_GPU_LJFLUID_BASE_GLUE_HPP
#define MDSIM_GPU_LJFLUID_BASE_GLUE_HPP

#include <cuda_wrapper.hpp>

namespace mdsim { namespace gpu { namespace ljfluid
{

extern cuda::symbol<unsigned int> npart;
extern cuda::symbol<float> box;
extern cuda::symbol<float> timestep;
extern cuda::symbol<float> r_cut;
extern cuda::symbol<float> rr_cut;
extern cuda::symbol<float> en_cut;

#ifdef USE_POTENTIAL_SMOOTHING
extern cuda::symbol<float> rri_smooth;
extern cuda::function <void (float3*, const float2)> sample_smooth_function;
#endif

#ifdef DIM_3D
extern cuda::function<void (float4*, float4*, float4*, float4 const*)> inteq;
extern cuda::function<void (float4*, unsigned int)> lattice;
extern cuda::function<void (float4*, unsigned int)> lattice_simple;
#else /* DIM_3D */
extern cuda::function<void (float2*, float2*, float2*, float2 const*)> inteq;
extern cuda::function<void (float2*, unsigned int)> lattice;
extern cuda::function<void (float2*, unsigned int)> lattice_simple;
#endif /* DIM_3D */
extern cuda::function<void (float const*, float2*)> potential_energy_sum;

}}} // namespace mdsim::gpu::ljfluid

#endif /* ! MDSIM_GPU_LJFLUID_BASE_GLUE_HPP */
