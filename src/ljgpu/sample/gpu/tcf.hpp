/* Time correlation functions for CUDA
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

#ifndef LJGPU_ALGORITHM_GPU_TCF_HPP
#define LJGPU_ALGORITHM_GPU_TCF_HPP

#include <cuda_wrapper.hpp>
#include <ljgpu/math/gpu/dsfun.cuh>

namespace ljgpu { namespace gpu
{

struct tcf_base
{
    enum {
	BLOCKS = 16,
	THREADS = 512,
	WARP_SIZE = 32,
    };

    static cuda::symbol<unsigned int*> count;
    static cuda::symbol<dfloat*> mean;
    static cuda::symbol<dfloat*> variance;
};

template <int dimension>
struct tcf;

template <>
struct tcf<3> : public tcf_base
{
    static cuda::function<void (float4 const*, float4 const*, uint)>
       	mean_square_displacement;
    static cuda::function<void (float4 const*, float4 const*, uint)>
       	mean_quartic_displacement;
    static cuda::function<void (float4 const*, float4 const*, uint)>
       	velocity_autocorrelation;
};

template <>
struct tcf<2> : public tcf_base
{
    static cuda::function<void (float2 const*, float2 const*, uint)>
       	mean_square_displacement;
    static cuda::function<void (float2 const*, float2 const*, uint)>
       	mean_quartic_displacement;
    static cuda::function<void (float2 const*, float2 const*, uint)>
       	velocity_autocorrelation;
};

}} // namespace ljgpu::gpu

#endif /* ! LJGPU_ALGORITHM_GPU_TCF_HPP */
