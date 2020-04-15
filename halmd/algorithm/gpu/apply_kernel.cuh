/*
 * Copyright © 2011  Michael Kopp
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <boost/preprocessor/repetition/enum_params.hpp>

#include <halmd/algorithm/gpu/apply_kernel.hpp>
#include <halmd/utility/gpu/thread.cuh>

namespace halmd {
namespace algorithm {
namespace gpu {
namespace apply_kernel {

/**
 * parallel application of a functor to some data
 */
template <
    typename functor
  , typename input_type
  , typename coalesced_input_type
  , typename output_type
  , typename coalesced_output_type
>
__global__ void apply(
    coalesced_input_type const* g_in
  , coalesced_output_type* g_out
  , unsigned int n
)
{
    // coalesced read- and write access to global memory for suitable array types
    for (unsigned int i = GTID; i < n; i += GTDIM) {
        g_out[i] = transform<functor, input_type, output_type>(g_in[i]);
    }
}

} // namespace apply_kernel

// bind function to wrapper

//
// To avoid repeating template arguments and at the same time not
// define an ugly macro, we use BOOST_PP_ENUM_PARAMS to generate
// the template arguments. The meaningless names (T0, T1, …) will
// never show up in compile error messages, as the compiler uses
// the template argument names of the *declaration*.
//

template <BOOST_PP_ENUM_PARAMS(5, typename T)>
apply_wrapper<BOOST_PP_ENUM_PARAMS(5, T)> apply_wrapper<BOOST_PP_ENUM_PARAMS(5, T)>::kernel = {
    apply_kernel::apply<BOOST_PP_ENUM_PARAMS(5, T)>
};

} // namespace algorithm
} // namespace gpu
} // namespace halmd
