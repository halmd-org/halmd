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

#ifndef LJGPU_MDSIM_GPU_LJFLUID_NBR_HPP
#define LJGPU_MDSIM_GPU_LJFLUID_NBR_HPP

#include <cuda_wrapper.hpp>

namespace ljgpu { namespace gpu { namespace ljfluid_neighbour
{

enum {
    /** fixed number of placeholders per cell */
    CELL_SIZE = 32,
    /** virtual particle tag */
    VIRTUAL_PARTICLE = -1,
};

extern cuda::function<void (float2*, float2*, float2*, float2 const*),
		      void (float4*, float4*, float4*, float4 const*)> inteq;
extern cuda::function<void (float3*, const float2)> sample_smooth_function;
extern cuda::function<void (float4 const*, float4*, float4*, int const*, float*, float*),
		      void (float2 const*, float2*, float2*, int const*, float*, float*)> mdstep;
extern cuda::function<void (int const*, int*, float4*),
		      void (int const*, int*, float2*)> update_neighbours;
extern cuda::function<void (float4 const*, uint*),
		      void (float2 const*, uint*)> compute_cell;
extern cuda::function<void (const int*, float4*, float4*, float4*, int*),
		      void (const int*, float2*, float2*, float2*, int*)> order_particles;
extern cuda::function<void (int*)> init_tags;
extern cuda::function<void (uint const*, int const*, int const*, int*)> assign_cells;
extern cuda::function<void (uint*, int*)> find_cell_offset;
extern cuda::function<void (int*)> gen_index;

extern cuda::symbol<uint> npart;
extern cuda::symbol<float> box;
extern cuda::symbol<float> timestep;
extern cuda::symbol<float> r_cut;
extern cuda::symbol<float> rr_cut;
extern cuda::symbol<float> en_cut;
extern cuda::symbol<float> rri_smooth;

extern cuda::symbol<uint> ncell;
extern cuda::symbol<uint> nbl_size;
extern cuda::symbol<uint> nbl_stride;
extern cuda::symbol<float> r_cell;
extern cuda::symbol<float> rr_cell;

template <typename T>
struct texref
{
    static cuda::texture<T> r;
    static cuda::texture<T> R;
    static cuda::texture<T> v;
    static cuda::texture<int> tag;
};

}}} // namespace ljgpu::gpu::ljfluid_neighbour

#endif /* ! LJGPU_MDSIM_GPU_LJFLUID_BASE_HPP */
