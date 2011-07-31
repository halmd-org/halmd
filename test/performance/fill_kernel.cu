/*
 * Copyright © 2011  Felix Höfling
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

#include <test/performance/fill_kernel.hpp>

// thread ID within block
#define TID     threadIdx.x
// number of threads per block
#define TDIM    blockDim.x
// block ID within grid
#define BID     (blockIdx.y * gridDim.x + blockIdx.x)
// number of blocks within grid
#define BDIM    (gridDim.y * gridDim.x)
// thread ID within grid
#define GTID    (BID * TDIM + TID)
// number of threads per grid
#define GTDIM   (BDIM * TDIM)

/**
 * fill float array with a given value using a loop
 */
__global__ void fill_loop(float* g_a, unsigned int size, float value)
{
    const unsigned int threads = GTDIM;

    for (unsigned int i = GTID; i < size; i += threads)
    {
        g_a[i] = value;
    }
}

/**
 * fill float array with a given value using just a single if,
 * this kernel requires GTDIM ≥ size
 */
__global__ void fill_if(float* g_a, unsigned int size, float value)
{
    const unsigned int i = GTID;
    if (i < size) {
        g_a[i] = value;
    }
}

/**
 * fill float array with a given value without any conditions,
 * this kernel requires GTDIM = capacity(g_a)
 */
__global__ void fill_all(float* g_a, float value)
{
    g_a[GTID] = value;
}

cuda::function<void (float*, unsigned int, float)> fill_loop_kernel(&fill_loop);
cuda::function<void (float*, unsigned int, float)> fill_if_kernel(&fill_if);
cuda::function<void (float*, float)> fill_all_kernel(&fill_all);
