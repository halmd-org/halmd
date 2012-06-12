/*
 * Copyright © 2008-2012 Peter Colberg
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

#ifndef HALMD_ALGORITHM_GPU_REDUCE_KERNEL_CUH
#define HALMD_ALGORITHM_GPU_REDUCE_KERNEL_CUH

#include <boost/utility/enable_if.hpp>

#include <halmd/algorithm/gpu/reduce_kernel.hpp>
#include <halmd/utility/gpu/thread.cuh>

namespace halmd {
namespace detail {

template <unsigned int threads, typename accumulator_type>
inline __device__ typename boost::enable_if_c<threads == 1, void>::type
reduce(accumulator_type& acc, accumulator_type s_acc[])
{
    if (TID < threads) {
        acc(s_acc[TID + threads]);
    }
}

template <unsigned int threads, typename accumulator_type>
inline __device__ typename boost::disable_if_c<threads == 1, void>::type
reduce(accumulator_type& acc, accumulator_type s_acc[])
{
    if (TID < threads) {
        acc(s_acc[TID + threads]);
        s_acc[TID] = acc;
    }

    //
    // FIXME __syncthreads only if threads >= warpSize
    //
    // On Fermi devices (even with PTX < 2.0), this requires acc and s_acc
    // to be declared as volatile, which in turn requires the methods of
    // the accumulator and any used data type (fixed_vector, dsfloat) to
    // be declared volatile, in the case of member methods, and/or accept
    // function arguments with the volatile qualifier.
    //
    // Refer to Fermi Compatibility Guide, 1.3.3 Kernels.
    //
    __syncthreads();

    reduce<threads / 2>(acc, s_acc);
}

/**
 * Reduce accumulators of block threads.
 *
 * @param acc accumulator of block thread 0
 */
template <unsigned int threads, typename accumulator_type>
inline __device__ void reduce(accumulator_type& acc)
{
    // We need to avoid default initialization of the shared memory
    // array, since this increases execution time of the kernel.
    // Use a char array, and cast to type of reduction functor.
    __shared__ char s_storage[threads * sizeof(accumulator_type)];

    accumulator_type* const s_acc = reinterpret_cast<accumulator_type*>(s_storage);

    // reduced value for this thread
    s_acc[TID] = acc;
    __syncthreads();

    // compute reduced value for all threads in block
    reduce<threads / 2>(acc, s_acc);
}

/**
 * Compute block sums of input array using unary accumulator.
 *
 * @param first input iterator to first element
 * @param size number of elements
 * @param g_block_acc output block accumulators
 * @param input accumulator
 */
template <unsigned int threads, typename accumulator_type>
static __global__ void reduction(
    typename accumulator_type::iterator const first
  , typename std::iterator_traits<typename accumulator_type::iterator>::difference_type size
  , accumulator_type* g_block_acc
  , accumulator_type acc
)
{
    // load values from global device memory
    for (unsigned int i = GTID; i < size; i += GTDIM) {
        acc(first[i]);
    }
    // compute reduced value for all threads in block
    reduce<threads>(acc);

    if (TID < 1) {
        // store block reduced value in global memory
        g_block_acc[blockIdx.x] = acc;
    }
}

template <unsigned int threads, typename accumulator_type>
reduction_kernel<threads, accumulator_type> const reduction_kernel<threads, accumulator_type>::kernel = {
    reduction<threads>
};

} // namespace detail
} // namespace halmd

#endif /* ! HALMD_ALGORITHM_GPU_REDUCE_KERNEL_CUH */
