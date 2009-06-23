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

#ifndef LJGPU_MDSIM_GPU_VELOCITY_HPP
#define LJGPU_MDSIM_GPU_VELOCITY_HPP

#include <cuda_wrapper.hpp>
#include <ljgpu/math/gpu/dsfloat.cuh>
#include <ljgpu/math/gpu/dsvector.cuh>

namespace ljgpu { namespace gpu { namespace velocity
{

enum {
    BLOCKS = 32,
    THREADS = 256,
};

extern cuda::function<void(float4 const*, dsfloat*, uint, uint), void(float2 const*, dsfloat*, uint, uint)> sum;

}}} // namespace ljgpu::gpu::velocity

#endif /* ! LJGPU_MDSIM_GPU_VELOCITY_HPP */
