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

#include <ljgpu/ljfluid/gpu/base.cuh>
#include <ljgpu/ljfluid/gpu/ljfluid_square.hpp>

namespace ljgpu { namespace gpu
{

/**
 * MD simulation step
 */
template <typename T, typename TT, typename U>
__global__ void mdstep(U* g_r, U* g_v, U* g_f, float* g_en, float* g_virial)
{
    extern __shared__ T s_r[];

    // load particle associated with this thread
    T r = unpack(g_r[GTID]);
    T v = unpack(g_v[GTID]);

    // potential energy contribution
    float en = 0;
    // virial equation sum contribution
    float virial = 0;
    // force sum
    TT f = 0;

    // iterate over all blocks
    for (unsigned int k = 0; k < gridDim.x; k++) {
	// load positions of particles within block
	s_r[TID] = unpack(g_r[k * blockDim.x + TID]);
	__syncthreads();

	// iterate over all particles within block
	for (unsigned int j = 0; j < blockDim.x; j++) {
	    // skip placeholder particles
	    if (k * blockDim.x + j >= npart)
		continue;
	    // skip identical particle
	    if (blockIdx.x == k && TID == j)
		continue;

	    // compute Lennard-Jones force with particle
	    compute_force(r, s_r[j], f, en, virial);
	}
	__syncthreads();
    }

    // second leapfrog step of integration of equations of motion
    leapfrog_full_step(v, f.f0);

    // store particle associated with this thread
    g_v[GTID] = pack(v);
    g_f[GTID] = pack(f.f0);
    g_en[GTID] = en;
    g_virial[GTID] = virial;
}

/**
 * device constant wrappers
 */
cuda::symbol<uint> ljfluid_base<ljfluid_impl_gpu_square>::npart(gpu::npart);
cuda::symbol<float> ljfluid_base<ljfluid_impl_gpu_square>::box(gpu::box);
cuda::symbol<float> ljfluid_base<ljfluid_impl_gpu_square>::timestep(gpu::timestep);
cuda::symbol<float> ljfluid_base<ljfluid_impl_gpu_square>::r_cut(gpu::r_cut);
cuda::symbol<float> ljfluid_base<ljfluid_impl_gpu_square>::rr_cut(gpu::rr_cut);
cuda::symbol<float> ljfluid_base<ljfluid_impl_gpu_square>::en_cut(gpu::en_cut);
cuda::symbol<float> ljfluid_base<ljfluid_impl_gpu_square>::rri_smooth(gpu::rri_smooth);

/**
 * device function wrappers
 */
cuda::function<void (float3*, const float2)>
    ljfluid_base<ljfluid_impl_gpu_square>::sample_smooth_function(gpu::sample_smooth_function);

cuda::function<void (float4*, float4*, float4*, float4 const*)>
    ljfluid<ljfluid_impl_gpu_square<3> >::inteq(gpu::inteq<float3>);
cuda::function<void (float4*, float4*, float4*, float*, float*)>
    ljfluid<ljfluid_impl_gpu_square<3> >::mdstep(gpu::mdstep<float3, dfloat3>);

cuda::function<void (float2*, float2*, float2*, float2 const*)>
    ljfluid<ljfluid_impl_gpu_square<2> >::inteq(gpu::inteq<float2>);
cuda::function<void (float2*, float2*, float2*, float*, float*)>
    ljfluid<ljfluid_impl_gpu_square<2> >::mdstep(gpu::mdstep<float2, dfloat2>);

}} // namespace ljgpu::gpu
