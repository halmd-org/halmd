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

#include <halmd/algorithm/gpu/apply_kernel.cuh>
#include <halmd/numeric/blas/blas.hpp>

using namespace halmd;
using namespace halmd::algorithm::gpu;

template class apply_wrapper<
    square_               // transform
  , fixed_vector<float, 2>  // input_type
  , float2                  // coalesced_input_type
  , float                   // output_type
  , float                   // coalesced_output_type
>;
