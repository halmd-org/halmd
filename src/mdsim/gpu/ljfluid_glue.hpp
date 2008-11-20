/* Lennard-Jones fluid kernel
 *
 * Copyright (C) 2008  Peter Colberg
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

#ifndef MDSIM_GPU_LJFLUID_GLUE_HPP
#define MDSIM_GPU_LJFLUID_GLUE_HPP

#include <cuda_wrapper.hpp>

/** virtual particle tag */
#define VIRTUAL_PARTICLE	-1

#ifdef USE_CELL

#ifndef CELL_SIZE
/** fixed number of placeholders per cell */
#define CELL_SIZE 32
#endif

/** real particle tag */
#define REAL_PARTICLE(x)	(x)

#define IS_REAL_PARTICLE(x)	(x >= 0)

#endif /* USE_CELL */

namespace mdsim { namespace gpu { namespace ljfluid
{

#ifdef DIM_3D
extern cuda::function<void (float4*, float4*, float4*, float4 const*)> inteq;
# ifdef USE_CELL
extern cuda::function<void (float4 const*, float4*, float4*, int const*, float*, float*)> mdstep;
extern cuda::function<void (float4 const*, float*)> maximum_velocity;
extern cuda::function<void (float4 const*, uint*)> compute_cell;
extern cuda::function<void (const int*, float4*, float4*, float4*, int*)> order_particles;
extern cuda::texture<float4> r;
extern cuda::texture<float4> R;
extern cuda::texture<float4> v;
# else
extern cuda::function<void (float4*, float4*, float4*, float*, float*)> mdstep;
# endif
extern cuda::function<void (float4*, unsigned int)> lattice;
extern cuda::function<void (float4*, unsigned int)> lattice_simple;
#else /* DIM_3D */
extern cuda::function<void (float2*, float2*, float2*, float2 const*)> inteq;
# ifdef USE_CELL
extern cuda::function<void (float2 const*, float2*, float2*, int const*, float*, float*)> mdstep;
extern cuda::function<void (float2 const*, float*)> maximum_velocity;
extern cuda::function<void (float2 const*, uint*)> compute_cell;
extern cuda::function<void (const int*, float2*, float2*, float2*, int*)> order_particles;
extern cuda::texture<float2> r;
extern cuda::texture<float2> R;
extern cuda::texture<float2> v;
# else
extern cuda::function<void (float2*, float2*, float2*, float*, float*)> mdstep;
# endif
extern cuda::function<void (float2*, unsigned int)> lattice;
extern cuda::function<void (float2*, unsigned int)> lattice_simple;
#endif /* DIM_3D */

extern cuda::symbol<unsigned int> npart;
extern cuda::symbol<float> box;
extern cuda::symbol<float> timestep;
extern cuda::symbol<float> r_cut;
extern cuda::symbol<float> rr_cut;
extern cuda::symbol<float> en_cut;
extern cuda::function<void (int*)> init_tags;
extern cuda::function<void (float const*, float2*)> potential_energy_sum;

#ifdef USE_CELL
extern cuda::function<void (int const*, int*)> update_neighbours;
extern cuda::symbol<unsigned int> ncell;
extern cuda::symbol<unsigned int> nbl_size;
extern cuda::symbol<unsigned int> nbl_stride;
extern cuda::symbol<float> r_cell;
extern cuda::symbol<float> rr_cell;
extern cuda::texture<int> tag;
extern cuda::function<void (uint const*, int const*, int const*, int*)> assign_cells;
extern cuda::function<void (uint*, int*)> find_cell_offset;
extern cuda::function<void (int*)> gen_index;
#endif

#ifdef USE_POTENTIAL_SMOOTHING
extern cuda::symbol<float> rri_smooth;
extern cuda::function <void (float3*, const float2)> sample_smooth_function;
#endif

}}} // namespace mdsim::gpu::ljfluid

#endif /* ! MDSIM_GPU_LJFLUID_GLUE_HPP */
