/*
 * Copyright © 2011  Michael Kopp
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

#include <halmd/algorithm/gpu/apply_kernel.cuh>
#include <halmd/numeric/blas/fixed_vector.hpp>

using namespace halmd;
using namespace halmd::algorithm::gpu;

// explicit instantiation of apply<negate> kernels

template class apply_wrapper<
    negate_                     // transform_functor
  , fixed_vector<float, 2>      // input_type
  , float4                      // coalesced_input_type
  , fixed_vector<float, 2>      // output_type
  , float4                      // coalesced_output_type
>;

template class apply_wrapper<
    negate_                     // transform_functor
  , fixed_vector<float, 3>      // input_type
  , float4                      // coalesced_input_type
  , fixed_vector<float, 3>      // output_type
  , float4                      // coalesced_output_type
>;
