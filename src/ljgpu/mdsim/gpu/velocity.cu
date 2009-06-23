/* Parallel reduction kernel
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

#include <ljgpu/algorithm/gpu/base.cuh>
#include <ljgpu/algorithm/gpu/reduce.cuh>
#include <ljgpu/math/gpu/dsvector.cuh>
#include <ljgpu/mdsim/gpu/velocity.hpp>

namespace ljgpu { namespace cu { namespace velocity
{

enum { THREADS = gpu::velocity::THREADS };

/**
 * squared velocity sum
 */
template <int dimension, typename T>
__global__ void sum(T const* g_v, dsfloat* g_block_sum, uint n, uint offset)
{
    __shared__ dsfloat s_vv[THREADS];
    dsfloat vv = 0;

    // load values from global device memory
    for (uint i = GTID; i < n; i += GTDIM) {
#ifdef USE_VERLET_DSFUN
	vector<dsfloat, dimension> v(g_v[i], g_v[i + offset]);
#else
	vector<float, dimension> v = g_v[i];
#endif
	vv += v * v;
    }
    // reduced value for this thread
    s_vv[TID] = vv;
    __syncthreads();

    // compute reduced value for all threads in block
    reduce<THREADS / 2, sum_>(vv, s_vv);

    if (TID < 1) {
	// store block reduced value in global memory
	g_block_sum[blockIdx.x] = vv;
    }
}

}}} // namespace ljgpu::cu::velocity

namespace ljgpu { namespace gpu
{

/**
 * device function wrappers
 */
cuda::function<void(float4 const*, dsfloat*, uint, uint), void(float2 const*, dsfloat*, uint, uint)>
    velocity::sum(cu::velocity::sum<3>, cu::velocity::sum<2>);

}} // namespace ljgpu::gpu
