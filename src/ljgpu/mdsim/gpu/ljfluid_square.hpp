/* Lennard-Jones fluid kernel
 *
 * Copyright © 2008-2009  Peter Colberg
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

#ifndef LJGPU_MDSIM_GPU_LJFLUID_SQUARE_HPP
#define LJGPU_MDSIM_GPU_LJFLUID_SQUARE_HPP

#include <cuda_wrapper.hpp>

namespace ljgpu { namespace gpu { namespace ljfluid_square
{

extern cuda::function<void (float2*, float2*, float2*, float2 const*),
		      void (float4*, float4*, float4*, float4 const*)> inteq;
extern cuda::function<void (float2*, uint),
		      void (float4*, uint)> lattice;
extern cuda::function<void (float2*, uint),
		      void (float4*, uint)> lattice_simple;
extern cuda::function<void (float const* g_en, float2* g_en_sum)> potential_energy_sum;
extern cuda::function<void (float3*, const float2)> sample_smooth_function;
extern cuda::function<void (float2*, float2*, float2*, float*, float*),
		      void (float4*, float4*, float4*, float*, float*)> mdstep;

extern cuda::symbol<uint> npart;
extern cuda::symbol<float> box;
extern cuda::symbol<float> timestep;
extern cuda::symbol<float> r_cut;
extern cuda::symbol<float> rr_cut;
extern cuda::symbol<float> en_cut;
extern cuda::symbol<float> rri_smooth;

}}} // namespace ljgpu::gpu::ljfluid_square

#endif /* ! LJGPU_MDSIM_GPU_LJFLUID_SQUARE_HPP */
