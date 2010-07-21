/*
 * Copyright © 2010  Peter Colberg
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

#include <halmd/algorithm/gpu/reduce_kernel.cuh>
#include <halmd/numeric/gpu/blas/dsfloat.cuh>
#include <halmd/numeric/gpu/blas/vector.cuh>

using namespace halmd::algorithm::gpu;
using namespace halmd::numeric::gpu::blas;

namespace halmd
{
namespace mdsim { namespace gpu
{
namespace thermodynamics_kernel
{

} // namespace thermodynamics_kernel

}} // namespace mdsim::gpu

template class reduce_wrapper<
    sum_                // reduce_transform
  , vector<float, 3>    // input_type
  , float4              // coalesced_input_type
  , dsfloat             // output_type
  , dsfloat             // coalesced_output_type
  , square_             // input_transform
>;

template class reduce_wrapper<
    sum_                // reduce_transform
  , vector<float, 2>    // input_type
  , float4              // coalesced_input_type
  , dsfloat             // output_type
  , dsfloat             // coalesced_output_type
  , square_             // input_transform
>;

template class reduce_wrapper<
    sum_                // reduce_transform
  , vector<float, 3>    // input_type
  , float4              // coalesced_input_type
  , vector<dsfloat, 3>  // output_type
  , vector<dsfloat, 3>  // coalesced_output_type
>;

template class reduce_wrapper<
    sum_                // reduce_transform
  , vector<float, 2>    // input_type
  , float4              // coalesced_input_type
  , vector<dsfloat, 2>  // output_type
  , vector<dsfloat, 2>  // coalesced_output_type
>;

template class reduce_wrapper<
    sum_                // reduce_transform
  , float               // input_type
  , float               // coalesced_input_type
  , dsfloat             // output_type
  , dsfloat             // coalesced_output_type
>;

} // namespace halmd
