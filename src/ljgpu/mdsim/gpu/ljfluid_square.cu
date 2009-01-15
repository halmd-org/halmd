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

#include <ljgpu/mdsim/gpu/base.cuh>
#include <ljgpu/mdsim/gpu/ljfluid_square.hpp>

namespace ljgpu { namespace cu { namespace ljfluid
{

/**
 * MD simulation step
 */
template <typename T, typename TT, ensemble_type ensemble, bool smooth, typename U>
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
	    compute_force<smooth>(r, s_r[j], f, en, virial);
	}
	__syncthreads();
    }

    // second leapfrog step of integration of equations of motion
    leapfrog_full_step(v, f.f0);
    // random collisions with heat bath
    if (ensemble == NVT) {
	anderson_thermostat(v);
    }

    // store particle associated with this thread
    g_v[GTID] = pack(v);
    g_f[GTID] = pack(f.f0);
    g_en[GTID] = en;
    g_virial[GTID] = virial;
}

}}} // namespace ljgpu::gpu::ljfluid

namespace ljgpu { namespace gpu
{

typedef ljfluid_base<ljfluid_impl_gpu_square> _Base;
typedef ljfluid<ljfluid_impl_gpu_square<3> > _3D;
typedef ljfluid<ljfluid_impl_gpu_square<2> > _2D;

/**
 * device function wrappers
 */
cuda::function<void (float4*, float4*, float4*, float*, float*)>
    _3D::mdstep(cu::ljfluid::mdstep<float3, dfloat3, cu::ljfluid::NVE, false>);
cuda::function<void (float4*, float4*, float4*, float*, float*)>
    _3D::mdstep_nvt(cu::ljfluid::mdstep<float3, dfloat3, cu::ljfluid::NVT, false>);
cuda::function<void (float4*, float4*, float4*, float*, float*)>
    _3D::mdstep_smooth(cu::ljfluid::mdstep<float3, dfloat3, cu::ljfluid::NVE, true>);
cuda::function<void (float4*, float4*, float4*, float*, float*)>
    _3D::mdstep_smooth_nvt(cu::ljfluid::mdstep<float3, dfloat3, cu::ljfluid::NVT, true>);

cuda::function<void (float2*, float2*, float2*, float*, float*)>
    _2D::mdstep(cu::ljfluid::mdstep<float2, dfloat2, cu::ljfluid::NVE, false>);
cuda::function<void (float2*, float2*, float2*, float*, float*)>
    _2D::mdstep_nvt(cu::ljfluid::mdstep<float2, dfloat2, cu::ljfluid::NVT, false>);
cuda::function<void (float2*, float2*, float2*, float*, float*)>
    _2D::mdstep_smooth(cu::ljfluid::mdstep<float2, dfloat2, cu::ljfluid::NVE, true>);
cuda::function<void (float2*, float2*, float2*, float*, float*)>
    _2D::mdstep_smooth_nvt(cu::ljfluid::mdstep<float2, dfloat2, cu::ljfluid::NVT, true>);

}} // namespace ljgpu::gpu
