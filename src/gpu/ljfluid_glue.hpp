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

#include <cuda_wrapper/cuda_wrapper.hpp>

namespace mdsim
{

/**
 * MD simulation state in global device memory
 */
template <typename T>
struct mdstep_param
{
    // particle coordinates
    T* r;
#ifndef USE_LEAPFROG
    // previous particle coordinates for Verlet algorithm
    T* rm;
#endif
    // particle velocities
    T* v;
    // particle forces
    T* f;
    // potential energies for each particle
    float* en;
    // virial equation sums for each particle
    float* virial;
};

} //namespace mdsim


namespace mdsim { namespace gpu { namespace ljfluid
{

#ifdef DIM_3D
extern cuda::function<void (mdstep_param<float3>)> mdstep;
extern cuda::function<void (float3*, float3, unsigned int)> init_lattice;
extern cuda::function<void (mdstep_param<float3>, float, ushort3*)> init_vel;
extern cuda::function<void (float3*)> init_forces;
#else
extern cuda::function<void (mdstep_param<float2>)> mdstep;
extern cuda::function<void (float2*, float2, unsigned int)> init_lattice;
extern cuda::function<void (mdstep_param<float2>, float, ushort3*)> init_vel;
extern cuda::function<void (float2*)> init_forces;
#endif

extern cuda::symbol<float> box;
extern cuda::symbol<float> timestep;
extern cuda::symbol<float> rr_cut;
extern cuda::symbol<float> en_cut;

extern cuda::symbol<uint3> a;
extern cuda::symbol<uint3> c;

}}} // namespace mdsim::gpu::ljfluid

#endif /* ! MDSIM_GPU_LJFLUID_GLUE_HPP */
