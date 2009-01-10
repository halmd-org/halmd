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

#ifndef LJGPU_LJFLUID_GPU_LJFLUID_SQUARE_HPP
#define LJGPU_LJFLUID_GPU_LJFLUID_SQUARE_HPP

#include <cuda_wrapper.hpp>
#include <ljgpu/ljfluid/impl.hpp>

namespace ljgpu { namespace gpu
{

template <template <int> class ljfluid_impl>
struct ljfluid_base;

template <>
struct ljfluid_base<ljfluid_impl_gpu_square>
{
    static cuda::symbol<uint> npart;
    static cuda::symbol<float> box;
    static cuda::symbol<float> timestep;
    static cuda::symbol<float> r_cut;
    static cuda::symbol<float> rr_cut;
    static cuda::symbol<float> en_cut;
    static cuda::symbol<float> rri_smooth;

    static cuda::function<void (float3*, const float2)> sample_smooth_function;
};

template <typename ljfluid_impl>
struct ljfluid;

template <>
struct ljfluid<ljgpu::ljfluid_impl_gpu_square<3> >
    : public ljfluid_base<ljfluid_impl_gpu_square>
{
    static cuda::function<void (float4*, float4*, float4*, float4 const*)> inteq;
    static cuda::function<void (float4*, float4*, float4*, float*, float*)> mdstep;
};

template <>
struct ljfluid<ljgpu::ljfluid_impl_gpu_square<2> >
    : public ljfluid_base<ljfluid_impl_gpu_square>
{
    static cuda::function<void (float2*, float2*, float2*, float2 const*)> inteq;
    static cuda::function<void (float2*, float2*, float2*, float*, float*)> mdstep;
};

}} // namespace ljgpu::gpu

#endif /* ! LJGPU_LJFLUID_GPU_LJFLUID_SQUARE_HPP */
