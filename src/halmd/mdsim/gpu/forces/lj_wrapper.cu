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

#include <halmd/mdsim/gpu/forces/lj_kernel.cu>
#include <halmd/mdsim/gpu/forces/lj_wrapper.cuh>

namespace halmd { namespace mdsim { namespace gpu { namespace forces
{

cuda::texture<float4>
  lj_wrapper<3>::r = lj_kernel::dim_<3>::r;
cuda::symbol<float3>
  lj_wrapper<3>::length = lj_kernel::dim_<3>::box_length;
cuda::symbol<unsigned int>
  lj_wrapper<3>::neighbor_size = lj_kernel::neighbor_size_;
cuda::symbol<unsigned int>
  lj_wrapper<3>::neighbor_stride = lj_kernel::neighbor_stride_;
cuda::texture<float4>
  lj_wrapper<3>::ljparam = lj_kernel::ljparam_;
cuda::function<void (float4*, unsigned int*, float*, float4*)>
  lj_wrapper<3>::compute = lj_kernel::compute<vector<float, 3> >;

cuda::texture<float4>
  lj_wrapper<2>::r = lj_kernel::dim_<2>::r;
cuda::symbol<float2>
  lj_wrapper<2>::length = lj_kernel::dim_<2>::box_length;
cuda::symbol<unsigned int>
  lj_wrapper<2>::neighbor_size = lj_kernel::neighbor_size_;
cuda::symbol<unsigned int>
  lj_wrapper<2>::neighbor_stride = lj_kernel::neighbor_stride_;
cuda::texture<float4>
  lj_wrapper<2>::ljparam = lj_kernel::ljparam_;
cuda::function<void (float2*, unsigned int*, float*, float2*)>
  lj_wrapper<2>::compute = lj_kernel::compute<vector<float, 2> >;

}}}} // namespace halmd::mdsim::gpu::forces
