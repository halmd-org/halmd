/* Maxwell-Boltzmann distribution at accurate temperature
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

#ifndef LJGPU_MDSIM_GPU_BOLTZMANN_HPP
#define LJGPU_MDSIM_GPU_BOLTZMANN_HPP

#include <cuda_wrapper.hpp>
#include <ljgpu/math/gpu/dsfloat.cuh>
#include <ljgpu/rng/gpu/uint48.cuh>

namespace ljgpu { namespace gpu
{

template <int dimension = 0>
struct boltzmann;

template <>
struct boltzmann<>
{
    enum { BLOCKS = 32 };
    enum { THREADS = 32 << DEVICE_SCALE };

    struct rand48
    {
	static cuda::symbol<uint48> a;
	static cuda::symbol<uint48> c;
	static cuda::symbol<ushort3*> state;
    };
};

template <>
struct boltzmann<3> : boltzmann<>
{
    static cuda::function<void (float4*, uint, uint, float, float4*)> gaussian;
    static cuda::function<void (float4*, uint, uint, float4 const*, dsfloat*)> shift_velocity;
    static cuda::function<void (float4*, uint, uint, dsfloat const*, dsfloat)> scale_velocity;
};

template <>
struct boltzmann<2> : boltzmann<>
{
    static cuda::function<void (float2*, uint, uint, float, float2*)> gaussian;
    static cuda::function<void (float2*, uint, uint, float2 const*, dsfloat*)> shift_velocity;
    static cuda::function<void (float2*, uint, uint, dsfloat const*, dsfloat)> scale_velocity;
};

}} // namespace ljgpu::gpu

#endif /* ! LJGPU_MDSIM_GPU_BOLTZMANN_HPP */
