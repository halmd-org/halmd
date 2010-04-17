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

#include <halmd/mdsim/gpu/integrator/verlet_kernel.cu>
#include <halmd/mdsim/gpu/integrator/verlet_wrapper.cuh>

namespace halmd { namespace mdsim { namespace gpu { namespace integrator
{

cuda::symbol<float>
  verlet_wrapper<3>::timestep = verlet_kernel::timestep_;
cuda::symbol<float3>
  verlet_wrapper<3>::length = verlet_kernel::box<3>::length;
cuda::function <void (float4*, float4*, float4*, float4 const*)>
  verlet_wrapper<3>::integrate = verlet_kernel::_integrate<vector<dsfloat, 3>, vector<float, 3> >;
cuda::function <void (float4*, float4 const*)>
  verlet_wrapper<3>::finalize = verlet_kernel::_finalize<vector<dsfloat, 3>, vector<float, 3> >;

cuda::symbol<float>
  verlet_wrapper<2>::timestep = verlet_kernel::timestep_;
cuda::symbol<float2>
  verlet_wrapper<2>::length = verlet_kernel::box<2>::length;
cuda::function <void (float4*, float2*, float4*, float2 const*)>
  verlet_wrapper<2>::integrate = verlet_kernel::_integrate<vector<dsfloat, 2>, vector<float, 2> >;
cuda::function <void (float4*, float2 const*)>
  verlet_wrapper<2>::finalize = verlet_kernel::_finalize<vector<dsfloat, 2>, vector<float, 2> >;

}}}} // namespace halmd::mdsim::gpu::integrator
