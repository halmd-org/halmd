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

#include <halmd/mdsim/gpu/sample/trajectory_kernel.cu>
#include <halmd/mdsim/gpu/sample/trajectory_wrapper.cuh>

namespace halmd { namespace mdsim { namespace gpu { namespace sample
{

cuda::texture<float4>
  trajectory_wrapper<3>::r = trajectory_kernel::dim_<3>::r;
cuda::texture<float4>
  trajectory_wrapper<3>::image = trajectory_kernel::dim_<3>::image;
cuda::texture<float4>
  trajectory_wrapper<3>::v = trajectory_kernel::dim_<3>::v;
cuda::symbol<float3>
  trajectory_wrapper<3>::box_length = trajectory_kernel::dim_<3>::box_length;
cuda::function<void (unsigned int const*, float4*, float4*)>
  trajectory_wrapper<3>::sample = trajectory_kernel::sample<vector<float, 3> >;

cuda::texture<float4>
  trajectory_wrapper<2>::r = trajectory_kernel::dim_<2>::r;
cuda::texture<float2>
  trajectory_wrapper<2>::image = trajectory_kernel::dim_<2>::image;
cuda::texture<float4>
  trajectory_wrapper<2>::v = trajectory_kernel::dim_<2>::v;
cuda::symbol<float2>
  trajectory_wrapper<2>::box_length = trajectory_kernel::dim_<2>::box_length;
cuda::function<void (unsigned int const*, float2*, float2*)>
  trajectory_wrapper<2>::sample = trajectory_kernel::sample<vector<float, 2> >;

}}}} // namespace halmd::mdsim::gpu::sample
