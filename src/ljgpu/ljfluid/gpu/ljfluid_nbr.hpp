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

#ifndef LJGPU_LJFLUID_GPU_LJFLUID_NBR_HPP
#define LJGPU_LJFLUID_GPU_LJFLUID_NBR_HPP

#include <cuda_wrapper.hpp>
#include <ljgpu/ljfluid/impl.hpp>

namespace ljgpu { namespace gpu
{

template <template <int> class ljfluid_impl>
struct ljfluid_base;

template <>
struct ljfluid_base<ljfluid_impl_gpu_neighbour>
{
    enum {
	/** fixed number of placeholders per cell */
	CELL_SIZE = 32,
	/** virtual particle tag */
	VIRTUAL_PARTICLE = -1,
    };

    static cuda::symbol<uint> npart;
    static cuda::symbol<float> box;
    static cuda::symbol<float> timestep;
    static cuda::symbol<float> r_cut;
    static cuda::symbol<float> rr_cut;
    static cuda::symbol<float> en_cut;
    static cuda::symbol<float> rri_smooth;

    static cuda::symbol<uint> ncell;
    static cuda::symbol<uint> nbl_size;
    static cuda::symbol<uint> nbl_stride;
    static cuda::symbol<float> r_cell;
    static cuda::symbol<float> rr_cell;

    static cuda::texture<int> tag;

    static cuda::function<void (float3*, const float2)> sample_smooth_function;
    static cuda::function<void (int*)> init_tags;
    static cuda::function<void (uint const*, int const*, int const*, int*)> assign_cells;
    static cuda::function<void (uint*, int*)> find_cell_offset;
    static cuda::function<void (int*)> gen_index;
};

template <typename ljfluid_impl>
struct ljfluid;

template <>
struct ljfluid<ljgpu::ljfluid_impl_gpu_neighbour<3> >
    : public ljfluid_base<ljfluid_impl_gpu_neighbour>
{
    static cuda::texture<float4> r;
    static cuda::texture<float4> R;
    static cuda::texture<float4> v;

    static cuda::function<void (float4*, float4*, float4*, float4 const*)> inteq;
    static cuda::function<void (float4 const*, float4*, float4*, int const*, float*, float*)> mdstep;
    static cuda::function<void (int const*, int*, float4*)> update_neighbours;
    static cuda::function<void (float4 const*, uint*)> compute_cell;
    static cuda::function<void (const int*, float4*, float4*, float4*, int*)> order_particles;
};

template <>
struct ljfluid<ljgpu::ljfluid_impl_gpu_neighbour<2> >
    : public ljfluid_base<ljfluid_impl_gpu_neighbour>
{
    static cuda::texture<float2> r;
    static cuda::texture<float2> R;
    static cuda::texture<float2> v;

    static cuda::function<void (float2*, float2*, float2*, float2 const*)> inteq;
    static cuda::function<void (float2 const*, float2*, float2*, int const*, float*, float*)> mdstep;
    static cuda::function<void (int const*, int*, float2*)> update_neighbours;
    static cuda::function<void (float2 const*, uint*)> compute_cell;
    static cuda::function<void (const int*, float2*, float2*, float2*, int*)> order_particles;
};

}} // namespace ljgpu::gpu

#endif /* ! LJGPU_LJFLUID_GPU_LJFLUID_BASE_HPP */
