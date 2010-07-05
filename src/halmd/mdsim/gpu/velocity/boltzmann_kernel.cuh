/* Maxwell-Boltzmann distribution at accurate temperature
 *
 * Copyright © 2008-2010  Peter Colberg and Felix Höfling
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
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

#ifndef HALMD_MDSIM_GPU_VELOCITY_BOLTZMANN_KERNEL_CUH
#define HALMD_MDSIM_GPU_VELOCITY_BOLTZMANN_KERNEL_CUH

#include <cuda_wrapper.hpp>
#include <halmd/numeric/gpu/blas/dsfloat.cuh>
#include <halmd/numeric/gpu/uint48.cuh>

namespace halmd
{
namespace mdsim { namespace gpu { namespace velocity
{

using numeric::gpu::uint48;
using numeric::gpu::blas::dsfloat;

template <int dimension = 0>
struct boltzmann_wrapper;

template <>
struct boltzmann_wrapper<>
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
struct boltzmann_wrapper<3> : boltzmann_wrapper<>
{
    static cuda::function<void (float4*, uint, uint, float, float4*)> gaussian;
    static cuda::function<void (float4*, uint, uint, float4 const*, dsfloat*)> shift_velocity;
    static cuda::function<void (float4*, uint, uint, dsfloat const*, dsfloat)> scale_velocity;
};

template <>
struct boltzmann_wrapper<2> : boltzmann_wrapper<>
{
    static cuda::function<void (float4*, uint, uint, float, float2*)> gaussian;
    static cuda::function<void (float4*, uint, uint, float2 const*, dsfloat*)> shift_velocity;
    static cuda::function<void (float4*, uint, uint, dsfloat const*, dsfloat)> scale_velocity;
};

}}} // namespace mdsim::gpu::velocity

} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_VELOCITY_BOLTZMANN_KERNEL_CUH */
