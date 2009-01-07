/* Bits and pieces used in GPU algorithms
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

#ifndef LJGPU_ALGORITHM_GPU_BASE_CUH
#define LJGPU_ALGORITHM_GPU_BASE_CUH

namespace ljgpu { namespace gpu
{

/**
 * swap arguments
 */
template <typename T>
__device__ __host__ void swap(T& a, T& b)
{
    T c = b;
    b = a;
    a = c;
}

/*
 * CUDA kernel helper macros
 */

// For thread configuration within a block, we only make use of 1
// dimension. While CUDA execution configuration allows (x, y, z)
// dimensions of (512, 512, 64) as a maximum, CUDA devices with
// compute capatibility 1.0 only support a total of 512 threads
// per block.
//

// thread ID within block
#define TID	threadIdx.x
// number of threads per block
#define TDIM	blockDim.x
// thread ID within grid
#define GTID	((blockIdx.y * gridDim.x + blockIdx.x) * TDIM + TID)
// number of threads per grid
#define GTDIM	(gridDim.y * gridDim.x * TDIM)


/*
 * CUDA kernel debugging
 */

#ifdef __DEVICE_EMULATION__
#include <stdio.h>
#include <assert.h>
#endif

}} // namespace ljgpu::gpu

#endif /* ! LJGPU_ALGORITHM_GPU_BASE_CUH */
