/*
 * Copyright © 2008-2010  Peter Colberg
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

#ifndef HALMD_SAMPLE_GPU_VACF_FILTER_HPP
#define HALMD_SAMPLE_GPU_VACF_FILTER_HPP

#include <cuda_wrapper.hpp>
#include <halmd/math/gpu/dsfloat.cuh>
#include <halmd/rng/gpu/uint48.cuh>

namespace halmd { namespace gpu
{

template <int dimension = 0>
struct vacf_filter;

template <>
struct vacf_filter<>
{
    enum { BLOCKS = 32 };
    enum { THREADS = 32 << DEVICE_SCALE };
};

template <>
struct vacf_filter<3> : vacf_filter<>
{
    static cuda::function<void (float4 const*, float4 const*, float4 const*, float4 const*, uint, float4*, float*)> accumulate;
    static cuda::function<void (float4 const*, float4 const*, uint, float const*, unsigned int*)> assign_index_from_msd;
};

template <>
struct vacf_filter<2> : vacf_filter<>
{
    static cuda::function<void (float2 const*, float2 const*, float2 const*, float2 const*, uint, float4*, float*)> accumulate;
    static cuda::function<void (float2 const*, float2 const*, uint, float const*, unsigned int*)> assign_index_from_msd;
};

}} // namespace halmd::gpu

#endif /* ! HALMD_SAMPLE_GPU_VACF_FILTER_HPP */
