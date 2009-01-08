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

#ifndef LJGPU_ALGORITHM_GPU_REDUCE_HPP
#define LJGPU_ALGORITHM_GPU_REDUCE_HPP

#include <cuda_wrapper.hpp>
#include <ljgpu/math/gpu/dsfun.cuh>

namespace ljgpu { namespace gpu { namespace reduce
{

enum {
    BLOCKS = 16,
    THREADS = 512,
};

extern cuda::function<void(float const*, dfloat*, uint)> sum;
extern cuda::function<void(float4 const*, float*, uint),
		      void(float2 const*, float*, uint)> max;

}}} // namespace ljgpu::gpu::reduce

#endif /* ! LJGPU_ALGORITHM_GPU_REDUCE_HPP */
