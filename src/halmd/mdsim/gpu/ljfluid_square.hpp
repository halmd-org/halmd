/* Lennard-Jones fluid kernel
 *
 * Copyright © 2008-2010  Peter Colberg
 *                        Felix Höfling
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

#ifndef HALMD_MDSIM_GPU_LJFLUID_SQUARE_HPP
#define HALMD_MDSIM_GPU_LJFLUID_SQUARE_HPP

#include <cuda_wrapper.hpp>
#include <halmd/mdsim/impl.hpp>
#include <halmd/mdsim/gpu/base.hpp>
#include <halmd/rng/gpu/uint48.cuh>

namespace halmd { namespace gpu
{

template <>
struct ljfluid_base<ljfluid_impl_gpu_square>
: public ljfluid_base<ljfluid_impl_gpu_base>
{
};

template <>
struct ljfluid<ljfluid_impl_gpu_square, 3>
: public ljfluid_base<ljfluid_impl_gpu_square>, public ljfluid<ljfluid_impl_gpu_base, 3>
{
    template <mixture_type, potential_type>
    struct variant
    {
        static cuda::function<void (float4 const*, float4*, float4*, float*, float4*, float4*)> mdstep;
    };

    static cuda::function<void (float4 const*, float4 const*, float4 const*, float4*, float4*)> sample;
    static cuda::function<void (float4*, float4*, float4*, float4 const*, float4*)> inteq;
};

template <>
struct ljfluid<ljfluid_impl_gpu_square, 2>
: public ljfluid_base<ljfluid_impl_gpu_square>, public ljfluid<ljfluid_impl_gpu_base, 2>
{
    template <mixture_type, potential_type>
    struct variant
    {
        static cuda::function<void (float4 const*, float2*, float2*, float*, float2*, float2*)> mdstep;
    };

    static cuda::function<void (float4 const*, float2 const*, float2 const*, float2*, float2*)> sample;
    static cuda::function<void (float4*, float2*, float2*, float2 const*, float2*)> inteq;
};

}} // namespace halmd::gpu

#endif /* ! HALMD_MDSIM_GPU_LJFLUID_SQUARE_HPP */
