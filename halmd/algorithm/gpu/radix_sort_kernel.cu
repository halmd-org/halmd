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

namespace halmd {
namespace algorithm {
namespace gpu {
namespace radix_sort_kernel {

/**
 * atomically add value to 32-bit word in shared memory
 */
template <unsigned int count>
__device__ void atomic_add(
    unsigned int const& i
  , unsigned int const& value
  , unsigned int& r
)
{
    extern __shared__ unsigned int s_bucket[];
    unsigned int const tid = threadIdx.x;

    // increment shared memory address within single thread of each half-warp
    if ((tid % HALF_WARP_SIZE) == (HALF_WARP_SIZE - count)) {
        r = s_bucket[i];
        s_bucket[i] = r + value;
    }
    __syncthreads();

    // recurse through all threads of each half-warp
    atomic_add<count - 1>(i, value, r);
}

template <>
__device__ void atomic_add<0>(
    unsigned int const&
  , unsigned int const&
  , unsigned int&
) {}

/**
 * returns 32-bit word in shared memory and atomically adds value
 */
__device__ unsigned int atomic_add(
    unsigned int const& i
  , unsigned int const& value
)
{
    unsigned int r;
    atomic_add<HALF_WARP_SIZE>(i, value, r);
    return r;
}

/**
 * compute partial radix counts for given radix shift
 */
__global__ void histogram_key(
    unsigned int const* g_input_key
  , unsigned int* g_bucket
  , unsigned int const count
  , unsigned int const shift
)
{
    //
    // Radix Sort for Vector Multiprocessors,
    // Marco Zagha and Guy E. Blelloch.
    // Supercomputing '91, November 1991. 
    //
    // http://www.cs.cmu.edu/~scandal/papers/cray-sort-supercomputing91.html
    //

    extern __shared__ unsigned int s_bucket[];

    unsigned int const tid = threadIdx.x;
    unsigned int const pid = threadIdx.x / HALF_WARP_SIZE;
    unsigned int const wid = threadIdx.x % HALF_WARP_SIZE;
    unsigned int const threads = blockDim.x;
    unsigned int const bid = blockIdx.x;
    unsigned int const blocks = gridDim.x;

    // set bucket counts to zero
    for (unsigned int i = 0; i < BUCKETS_PER_THREAD; ++i) {
        s_bucket[tid + i * threads] = 0;
    }
    __syncthreads();

    // number of partitions per block
    unsigned int const parts = threads / HALF_WARP_SIZE;
    // number of elements per partition, aligned to total thread count
    unsigned int const elems = ((count + blocks * threads - 1) / (blocks * threads)) * HALF_WARP_SIZE;

    for (unsigned int i = 0; i < elems; i += HALF_WARP_SIZE) {
        // global memory offset of sort key
        unsigned int const j = i + wid + (pid + parts * bid) * elems;
        // read sort key and transform according to radix shift
        unsigned int const radix = (j < count) ? (g_input_key[j] >> shift) & (BUCKET_SIZE - 1) : 0;

        // atomically increment bucket count
        atomic_add(pid + parts * radix, (j < count) ? 1 : 0);
    }

    // write radix counts to global memory
    for (unsigned int i = 0; i < BUCKETS_PER_THREAD; ++i) {
        // partition
        unsigned int const j = tid % parts;
        // bucket
        unsigned int const k = tid / parts + i * threads / parts;
        // write count to partition bucket in column major order
        g_bucket[j + (bid + k * blocks) * parts] = s_bucket[tid + i * threads];
    }
}

/**
 * permute array given radix counts prefix sums
 */
template <bool value>
__device__ void permute(
    unsigned int const* g_input_key
  , unsigned int* g_output_key
  , unsigned int const* g_bucket
  , unsigned int const count
  , unsigned int const shift
  , unsigned int const* g_input_value = 0
  , unsigned int* g_output_value = 0
)
{
    extern __shared__ unsigned int s_bucket[];

    unsigned int const tid = threadIdx.x;
    unsigned int const pid = threadIdx.x / HALF_WARP_SIZE;
    unsigned int const wid = threadIdx.x % HALF_WARP_SIZE;
    unsigned int const threads = blockDim.x;
    unsigned int const bid = blockIdx.x;
    unsigned int const blocks = gridDim.x;

    // number of partitions per block
    unsigned int const parts = threads / HALF_WARP_SIZE;
    // number of elements per partition, aligned to total thread count
    unsigned int const elems = ((count + blocks * threads - 1) / (blocks * threads)) * HALF_WARP_SIZE;

    // read radix counts from global memory
    for (unsigned int i = 0; i < BUCKETS_PER_THREAD; ++i) {
        // partition
        unsigned int const j = tid % parts;
        // bucket
        unsigned int const k = tid / parts + i * threads / parts;
        // read count from partition bucket in column major order
        s_bucket[tid + i * threads] = g_bucket[j + (bid + k * blocks) * parts];
    }
    __syncthreads();

    for (unsigned int i = 0; i < elems; i += HALF_WARP_SIZE) {
        // global memory offset of sort key
        unsigned int const j = i + wid + (pid + parts * bid) * elems;
        // read sort key from global memory
        unsigned int const key = (j < count) ? g_input_key[j] : 0;
        // transform sort key according to radix shift
        unsigned int const radix = (key >> shift) & (BUCKET_SIZE - 1);

        // atomically read and increment global radix offset
        unsigned int const l = atomic_add(pid + parts * radix, (j < count) ? 1 : 0);

        // scatter write permuted array element to global memory
        if (j < count) {
            // write sort key
            g_output_key[l] = key;
            // permute data array element
            if (value) {
                g_output_value[l] = g_input_value[j];
            }
        }
    }
}

__global__ void permute_key(
    unsigned int const* g_input_key
  , unsigned int* g_output_key
  , unsigned int const* g_bucket
  , unsigned int const count
  , unsigned int const shift
)
{
    permute<false>(
        g_input_key
      , g_output_key
      , g_bucket
      , count
      , shift
    );
}

__global__ void permute_key_value(
    unsigned int const* g_input_key
  , unsigned int* g_output_key
  , unsigned int const* g_bucket
  , unsigned int const count
  , unsigned int const shift
  , unsigned int const* g_input_value
  , unsigned int* g_output_value
)
{
    permute<true>(
        g_input_key
      , g_output_key
      , g_bucket
      , count
      , shift
      , g_input_value
      , g_output_value
    );
}

} // namespace radix_sort_kernel

/**
 * device function wrappers
 */
radix_sort_wrapper const radix_sort_wrapper::kernel = {
    radix_sort_kernel::histogram_key
  , radix_sort_kernel::permute_key
  , radix_sort_kernel::permute_key_value
};

template class scan_wrapper<unsigned int>;

} // namespace gpu
} // namespace algorithm
} // namespace halmd
