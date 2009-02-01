/* Parallel reduction
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

#ifndef LJGPU_ALGORITHM_GPU_REDUCE_CUH
#define LJGPU_ALGORITHM_GPU_REDUCE_CUH

#include <boost/mpl/int.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>
using namespace boost;

namespace ljgpu { namespace cu
{

enum { WARP_SIZE = 32 };

template <int threads, typename T>
__device__ typename enable_if<is_same<mpl::int_<threads>, mpl::int_<1> >, void>::type
reduce(T& sum, T s_sum[])
{
    if (threadIdx.x < threads) {
	sum += s_sum[threadIdx.x + threads];
    }
}

template <int threads, typename T>
__device__ typename disable_if<is_same<mpl::int_<threads>, mpl::int_<1> >, void>::type
reduce(T& sum, T s_sum[])
{
    if (threadIdx.x < threads) {
	sum += s_sum[threadIdx.x + threads];
	s_sum[threadIdx.x] = sum;
    }
    // no further syncs needed within execution warp of 32 threads
    if (threads >= WARP_SIZE) {
	__syncthreads();
    }

    reduce<threads / 2>(sum, s_sum);
}

}} // namespace ljgpu::cu

#endif /* LJGPU_ALGORITHM_GPU_REDUCE_CUH */
