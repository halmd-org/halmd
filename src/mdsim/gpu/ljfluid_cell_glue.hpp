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

#ifndef MDSIM_GPU_LJFLUID_CELL_GLUE_HPP
#define MDSIM_GPU_LJFLUID_CELL_GLUE_HPP

#include <cuda_wrapper.hpp>
#include "ljfluid_base_glue.hpp"

#ifndef CELL_SIZE
/** fixed number of placeholders per cell */
#define CELL_SIZE 32
#endif

/** real particle tag */
#define REAL_PARTICLE(x)	(x)
/** virtual particle tag */
#define VIRTUAL_PARTICLE	-1

#define IS_REAL_PARTICLE(x)	(x >= 0)

namespace mdsim { namespace gpu { namespace ljfluid_gpu_cell
{

extern cuda::symbol<unsigned int> ncell;

#ifdef DIM_3D
extern cuda::function<void (float4 const*, float4*, float4*, int const*, float*, float*)> mdstep;
extern cuda::function<void (float4 const*, float4*, int*)> assign_cells;
extern cuda::function<void (float4 const*, float4 const*, float4 const*, int const*, float4*, float4*, float4*, int*)> update_cells;
#else
extern cuda::function<void (float2 const*, float2*, float2*, int const*, float*, float*)> mdstep;
extern cuda::function<void (float2 const*, float2*, int*)> assign_cells;
extern cuda::function<void (float2 const*, float2 const*, float2 const*, int const*, float2*, float2*, float2*, int*)> update_cells;
#endif

}}} // namespace mdsim::gpu::ljfluid

#endif /* ! MDSIM_GPU_LJFLUID_CELL_GLUE_HPP */
