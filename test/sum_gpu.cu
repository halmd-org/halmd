/* Test alternatives for summing over a CUDA device vector
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

#include "sum_gpu.hpp"
using namespace cuda;

#ifdef __DEVICE_EMULATION__
#include "stdio.h"
#endif


namespace mdsim { namespace test {

/**
 * calculate sum of vector components for each block
 */
__global__ void blocksum(float* v, float* v_sum)
{
    extern __shared__ float vblock[];

    vblock[threadIdx.x] = v[blockDim.x * blockIdx.x + threadIdx.x];

    //
    // requires 512 threads per block
    //

    __syncthreads();
    if (0 == threadIdx.x % 2) vblock[threadIdx.x] += vblock[threadIdx.x + 1];
    __syncthreads();
    if (0 == threadIdx.x % 4) vblock[threadIdx.x] += vblock[threadIdx.x + 2];
    __syncthreads();
    if (0 == threadIdx.x % 8) vblock[threadIdx.x] += vblock[threadIdx.x + 4];
    __syncthreads();
    if (0 == threadIdx.x % 16) vblock[threadIdx.x] += vblock[threadIdx.x + 8];
    __syncthreads();
    if (0 == threadIdx.x % 32) vblock[threadIdx.x] += vblock[threadIdx.x + 16];
    __syncthreads();
    if (0 == threadIdx.x % 64) vblock[threadIdx.x] += vblock[threadIdx.x + 32];
    __syncthreads();
    if (0 == threadIdx.x % 128) vblock[threadIdx.x] += vblock[threadIdx.x + 64];
    __syncthreads();
    if (0 == threadIdx.x % 256) vblock[threadIdx.x] += vblock[threadIdx.x + 128];
    __syncthreads();
    if (0 == threadIdx.x % 512) vblock[threadIdx.x] += vblock[threadIdx.x + 256];
    __syncthreads();

    if (0 == threadIdx.x) {
	v_sum[blockIdx.x] = vblock[threadIdx.x];
    }
}

}} //namespace mdsim::test


namespace mdsim { namespace test { namespace gpu {

function<void (float*, float*)> blocksum(test::blocksum);

}}} // namespace mdsim::test::gpu
