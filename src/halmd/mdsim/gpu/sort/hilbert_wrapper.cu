/*
 * Copyright © 2008-2010  Peter Colberg
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

#include <halmd/mdsim/gpu/sort/hilbert_kernel.cu>
#include <halmd/mdsim/gpu/sort/hilbert_wrapper.cuh>

namespace halmd { namespace mdsim { namespace gpu
{

template <int dimension> cuda::symbol<float>
    hilbert_wrapper<dimension>::depth = hilbert_kernel::depth_;
template <int dimension> cuda::symbol<typename hilbert_wrapper<dimension>::vector_type>
    hilbert_wrapper<dimension>::box_length = hilbert_kernel::dim_<dimension>::box_length;

template <> cuda::function<void (float4 const*, unsigned int*)>
    hilbert_wrapper<3>::map = hilbert_kernel::map<vector<float, 3> >;
template <> cuda::function<void (float4 const*, unsigned int*)>
    hilbert_wrapper<2>::map = hilbert_kernel::map<vector<float, 2> >;

// explicit instantiation
template class hilbert_wrapper<3>;
template class hilbert_wrapper<2>;

}}} // namespace halmd::mdsim::gpu
