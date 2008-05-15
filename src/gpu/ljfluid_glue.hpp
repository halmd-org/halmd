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

#ifndef MDSIM_GPU_LJFLUID_GLUE_HPP
#define MDSIM_GPU_LJFLUID_GLUE_HPP

#include <cuda_wrapper.hpp>


namespace mdsim { namespace gpu { namespace ljfluid
{

#ifdef DIM_3D
extern cuda::function<void (float3*, float3*, float3*)> inteq;
extern cuda::function<void (float3*, float, ushort3*)> boltzmann;
extern cuda::function<void (float3*, float3*, float3*, float*, float*)> mdstep;
extern cuda::function<void (float3*)> lattice;
#else
extern cuda::function<void (float2*, float2*, float2*)> inteq;
extern cuda::function<void (float2*, float, ushort3*)> boltzmann;
extern cuda::function<void (float2*, float2*, float2*, float*, float*)> mdstep;
extern cuda::function<void (float2*)> lattice;
#endif

extern cuda::symbol<unsigned int> npart;
extern cuda::symbol<float> box;
extern cuda::symbol<float> timestep;
extern cuda::symbol<float> rr_cut;
extern cuda::symbol<float> en_cut;

extern cuda::symbol<uint3> a;
extern cuda::symbol<uint3> c;

}}} // namespace mdsim::gpu::ljfluid

#endif /* ! MDSIM_GPU_LJFLUID_GLUE_HPP */
