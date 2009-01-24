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
using namespace ljgpu::gpu;

namespace ljgpu { namespace cu { namespace ljfluid
{

/**
 * MD simulation step
 */
template <typename vector_type,
          mixture_type mixture,
	  potential_type potential,
	  ensemble_type ensemble,
	  typename T>
__global__ void mdstep(float4 const* g_r, T* g_v, T* g_f, float* g_en, float* g_virial)
{
    enum { dimension = vector_type::static_size };

    extern __shared__ unsigned int s_tag[];
    vector_type* const s_r = reinterpret_cast<vector_type*>(&s_tag[TDIM]);

    // load particle associated with this thread
    vector_type r, v;
    unsigned int tag;
    (r, tag) = g_r[GTID];
    v = g_v[GTID];
    // particle type in binary mixture
    int const a = (tag >= mpart[0]);

    // potential energy contribution
    float en = 0;
    // virial equation sum contribution
    float virial = 0;
    // force sum
    vector<dfloat, dimension> f = 0;

    // iterate over all blocks
    for (unsigned int k = 0; k < gridDim.x; k++) {
	// load positions of particles within block
	__syncthreads();
	(s_r[TID], s_tag[TID]) = g_r[k * blockDim.x + TID];
	__syncthreads();

	// iterate over all particles within block
	for (unsigned int j = 0; j < blockDim.x; j++) {
	    // skip placeholder particles
	    if (k * blockDim.x + j >= npart)
		continue;
	    // skip identical particle
	    if (blockIdx.x == k && TID == j)
		continue;

	    // particle type in binary mixture
	    int const b = (s_tag[j] >= mpart[0]);
	    // compute Lennard-Jones force with particle
	    compute_force<mixture, potential>(r, s_r[j], f, en, virial, a + b);
	}
    }

    // second leapfrog step of integration of equations of motion
    leapfrog_full_step(v, static_cast<vector_type>(f));

    // random collisions with heat bath
    if (ensemble == NVT) {
	anderson_thermostat(v);
    }

    // store particle associated with this thread
    g_v[GTID] = v;
    g_f[GTID] = static_cast<vector_type>(f);
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
cuda::function<void (float4 const*, float4*, float4*, float*, float*)>
    _3D::template variant<UNARY, C0POT, NVE>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 3>, UNARY, C0POT, NVE>);
cuda::function<void (float4 const*, float4*, float4*, float*, float*)>
    _3D::template variant<UNARY, C0POT, NVT>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 3>, UNARY, C0POT, NVT>);
cuda::function<void (float4 const*, float4*, float4*, float*, float*)>
    _3D::template variant<UNARY, C2POT, NVE>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 3>, UNARY, C2POT, NVE>);
cuda::function<void (float4 const*, float4*, float4*, float*, float*)>
    _3D::template variant<UNARY, C2POT, NVT>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 3>, UNARY, C2POT, NVT>);

cuda::function<void (float4 const*, float4*, float4*, float*, float*)>
    _3D::template variant<BINARY, C0POT, NVE>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 3>, BINARY, C0POT, NVE>);
cuda::function<void (float4 const*, float4*, float4*, float*, float*)>
    _3D::template variant<BINARY, C0POT, NVT>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 3>, BINARY, C0POT, NVT>);
cuda::function<void (float4 const*, float4*, float4*, float*, float*)>
    _3D::template variant<BINARY, C2POT, NVE>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 3>, BINARY, C2POT, NVE>);
cuda::function<void (float4 const*, float4*, float4*, float*, float*)>
    _3D::template variant<BINARY, C2POT, NVT>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 3>, BINARY, C2POT, NVT>);

cuda::function<void (float4 const*, float2*, float2*, float*, float*)>
    _2D::template variant<UNARY, C0POT, NVE>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 2>, UNARY, C0POT, NVE>);
cuda::function<void (float4 const*, float2*, float2*, float*, float*)>
    _2D::template variant<UNARY, C0POT, NVT>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 2>, UNARY, C0POT, NVT>);
cuda::function<void (float4 const*, float2*, float2*, float*, float*)>
    _2D::template variant<UNARY, C2POT, NVE>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 2>, UNARY, C2POT, NVE>);
cuda::function<void (float4 const*, float2*, float2*, float*, float*)>
    _2D::template variant<UNARY, C2POT, NVT>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 2>, UNARY, C2POT, NVT>);

cuda::function<void (float4 const*, float2*, float2*, float*, float*)>
    _2D::template variant<BINARY, C0POT, NVE>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 2>, BINARY, C0POT, NVE>);
cuda::function<void (float4 const*, float2*, float2*, float*, float*)>
    _2D::template variant<BINARY, C0POT, NVT>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 2>, BINARY, C0POT, NVT>);
cuda::function<void (float4 const*, float2*, float2*, float*, float*)>
    _2D::template variant<BINARY, C2POT, NVE>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 2>, BINARY, C2POT, NVE>);
cuda::function<void (float4 const*, float2*, float2*, float*, float*)>
    _2D::template variant<BINARY, C2POT, NVT>::mdstep(cu::ljfluid::mdstep<cu::vector<float, 2>, BINARY, C2POT, NVT>);

}} // namespace ljgpu::gpu
