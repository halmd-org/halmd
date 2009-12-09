/* Lennard-Jones fluid kernel
 *
 * Copyright © 2008-2009  Peter Colberg
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

#ifndef HALMD_MDSIM_GPU_LJFLUID_NBR_HPP
#define HALMD_MDSIM_GPU_LJFLUID_NBR_HPP

#include <cuda_wrapper.hpp>
#include <halmd/mdsim/impl.hpp>
#include <halmd/mdsim/gpu/base.hpp>
#include <halmd/rng/gpu/uint48.cuh>

namespace halmd { namespace gpu
{

template <>
struct ljfluid_base<ljfluid_impl_gpu_neighbour>
: public ljfluid_base<ljfluid_impl_gpu_base>
{
    static cuda::symbol<uint> ncell;
    static cuda::symbol<uint> nbl_size;
    static cuda::symbol<uint> nbl_stride;
    static cuda::symbol<float> rr_nbl;
    static cuda::symbol<unsigned int*> g_nbl;

    static cuda::function<void (unsigned int*, uint const*, unsigned int const*, unsigned int const*, unsigned int*)> assign_cells;
    static cuda::function<void (uint*, unsigned int*)> find_cell_offset;
    static cuda::function<void (unsigned int*)> gen_index;
};

template <>
struct ljfluid<halmd::ljfluid_impl_gpu_neighbour, 3>
: public ljfluid_base<ljfluid_impl_gpu_neighbour>, public ljfluid<ljfluid_impl_gpu_base, 3>
{
    static cuda::texture<float4> r;
    static cuda::texture<float4> R;
    static cuda::texture<float4> v;

    static cuda::function<void (unsigned int*, unsigned int const*)> update_neighbours;
    static cuda::function<void (float4 const*, uint*)> compute_cell;
    static cuda::function<void (unsigned int const*, float4*, float4*, float4*, unsigned int*)> order_particles;
    static cuda::function<void (uint const*, float4*)> order_velocities;
    static cuda::function<void (unsigned int const*, float4*, float4*)> sample;
    static cuda::function<void (float4*, float4*, float4*, float4*, float4 const*)> inteq;

    template <mixture_type, potential_type>
    struct variant
    {
        static cuda::function<void (float4 const*, float4*, float4*, float*, float4*)> mdstep;
    };
};

template <>
struct ljfluid<halmd::ljfluid_impl_gpu_neighbour, 2>
: public ljfluid_base<ljfluid_impl_gpu_neighbour>, public ljfluid<ljfluid_impl_gpu_base, 2>
{
    static cuda::texture<float4> r;
    static cuda::texture<float2> R;
    static cuda::texture<float2> v;

    static cuda::function<void (unsigned int*, unsigned int const*)> update_neighbours;
    static cuda::function<void (float4 const*, uint*)> compute_cell;
    static cuda::function<void (unsigned int const*, float4*, float2*, float2*, unsigned int*)> order_particles;
    static cuda::function<void (uint const*, float2*)> order_velocities;
    static cuda::function<void (unsigned int const*, float2*, float2*)> sample;
    static cuda::function<void (float4*, float2*, float2*, float2*, float2 const*)> inteq;

    template <mixture_type, potential_type>
    struct variant
    {
        static cuda::function<void (float4 const*, float2*, float2*, float*, float2*)> mdstep;
    };
};

}} // namespace halmd::gpu

#endif /* ! HALMD_MDSIM_GPU_LJFLUID_BASE_HPP */
