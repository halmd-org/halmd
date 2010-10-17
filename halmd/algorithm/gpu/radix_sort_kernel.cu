/*
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

#include <halmd/algorithm/gpu/radix_sort_kernel.hpp>
#include <halmd/algorithm/gpu/scan_kernel.cuh>

namespace halmd
{
namespace algorithm { namespace gpu
{
namespace radix_sort_kernel
{

/**
 * atomically add value to 32-bit word in shared memory
 */
template <uint count>
__device__ void atomic_add(uint const& i, uint const& value, uint& r)
{
    extern __shared__ uint s_bucket[];
    const uint tid = threadIdx.x;

    // increment shared memory address within single thread of each half-warp
    if ((tid % HALF_WARP_SIZE) == (HALF_WARP_SIZE - count)) {
        r = s_bucket[i];
        s_bucket[i] = r + value;
    }
    __syncthreads();

    // recurse through all threads of each half-warp
    atomic_add<count - 1>(i, value, r);
}

template <> __device__ void atomic_add<0>(uint const&, uint const&, uint&) {}

/**
 * returns 32-bit word in shared memory and atomically adds value
 */
__device__ uint atomic_add(uint const& i, uint const& value)
{
    uint r;
    atomic_add<HALF_WARP_SIZE>(i, value, r);
    return r;
}

/**
 * compute partial radix counts for given radix shift
 */
__global__ void histogram_keys(uint const* g_in, uint* g_bucket, const uint count, const uint shift)
{
    //
    // Radix Sort for Vector Multiprocessors,
    // Marco Zagha and Guy E. Blelloch.
    // Supercomputing '91, November 1991. 
    //
    // http://www.cs.cmu.edu/~scandal/papers/cray-sort-supercomputing91.html
    //

    extern __shared__ uint s_bucket[];

    const uint tid = threadIdx.x;
    const uint pid = threadIdx.x / HALF_WARP_SIZE;
    const uint wid = threadIdx.x % HALF_WARP_SIZE;
    const uint threads = blockDim.x;
    const uint bid = blockIdx.x;
    const uint blocks = gridDim.x;

    // set bucket counts to zero
    for (uint i = 0; i < BUCKETS_PER_THREAD; ++i) {
        s_bucket[tid + i * threads] = 0;
    }
    __syncthreads();

    // number of partitions per block
    const uint parts = threads / HALF_WARP_SIZE;
    // number of elements per partition, aligned to total thread count
    const uint elems = ((count + blocks * threads - 1) / (blocks * threads)) * HALF_WARP_SIZE;

    for (uint i = 0; i < elems; i += HALF_WARP_SIZE) {
        // global memory offset of sort key
        const uint j = i + wid + (pid + parts * bid) * elems;
        // read sort key and transform according to radix shift
        const uint radix = (j < count) ? (g_in[j] >> shift) & (BUCKET_SIZE - 1) : 0;

        // atomically increment bucket count
        atomic_add(pid + parts * radix, (j < count) ? 1 : 0);
    }

    // write radix counts to global memory
    for (uint i = 0; i < BUCKETS_PER_THREAD; ++i) {
        // partition
        const uint j = tid % parts;
        // bucket
        const uint k = tid / parts + i * threads / parts;
        // write count to partition bucket in column major order
        g_bucket[j + (bid + k * blocks) * parts] = s_bucket[tid + i * threads];
    }
}

/**
 * permute array given radix counts prefix sums
 */
template <typename T>
__global__ void permute(uint const* g_in, uint* g_out, T const* g_data_in, T* g_data_out, uint const* g_bucket, const uint count, const uint shift)
{
    extern __shared__ uint s_bucket[];

    const uint tid = threadIdx.x;
    const uint pid = threadIdx.x / HALF_WARP_SIZE;
    const uint wid = threadIdx.x % HALF_WARP_SIZE;
    const uint threads = blockDim.x;
    const uint bid = blockIdx.x;
    const uint blocks = gridDim.x;

    // number of partitions per block
    const uint parts = threads / HALF_WARP_SIZE;
    // number of elements per partition, aligned to total thread count
    const uint elems = ((count + blocks * threads - 1) / (blocks * threads)) * HALF_WARP_SIZE;

    // read radix counts from global memory
    for (uint i = 0; i < BUCKETS_PER_THREAD; ++i) {
        // partition
        const uint j = tid % parts;
        // bucket
        const uint k = tid / parts + i * threads / parts;
        // read count from partition bucket in column major order
        s_bucket[tid + i * threads] = g_bucket[j + (bid + k * blocks) * parts];
    }
    __syncthreads();

    for (uint i = 0; i < elems; i += HALF_WARP_SIZE) {
        // global memory offset of sort key
        const uint j = i + wid + (pid + parts * bid) * elems;
        // read sort key from global memory
        const uint key = (j < count) ? g_in[j] : 0;
        // transform sort key according to radix shift
        const uint radix = (key >> shift) & (BUCKET_SIZE - 1);

        // atomically read and increment global radix offset
        const uint l = atomic_add(pid + parts * radix, (j < count) ? 1 : 0);

        // scatter write permuted array element to global memory
        if (j < count) {
            // write sort key
            g_out[l] = key;
            // permute data array element
            g_data_out[l] = g_data_in[j];
        }
    }
}

} // namespace radix_sort_kernel

/**
 * device function wrappers
 */
template <typename T>
radix_sort_wrapper<T> const radix_sort_wrapper<T>::kernel = {
    radix_sort_kernel::histogram_keys
  , radix_sort_kernel::permute
};

template class radix_sort_wrapper<int>;
template class radix_sort_wrapper<unsigned int>;
template class radix_sort_wrapper<float4>;
template class radix_sort_wrapper<float2>;

}} // namespace algorithm::gpu

template class algorithm::gpu::scan_wrapper<unsigned int>;

} // namespace halmd
