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

#include <halmd/mdsim/gpu/neighbour_kernel.cu>
#include <halmd/mdsim/gpu/neighbour_wrapper.cuh>

namespace halmd { namespace mdsim { namespace gpu
{

cuda::texture<float>
  neighbour_wrapper<3>::rr_cut_skin = neighbour_kernel::rr_cut_skin_;
cuda::symbol<uint3>
  neighbour_wrapper<3>::ncell = neighbour_kernel::dim_<3>::ncell;
cuda::symbol<unsigned int>
  neighbour_wrapper<3>::neighbour_size = neighbour_kernel::neighbour_size_;
cuda::symbol<unsigned int>
  neighbour_wrapper<3>::neighbour_stride = neighbour_kernel::neighbour_stride_;
cuda::symbol<unsigned int>
  neighbour_wrapper<3>::nbox = neighbour_kernel::nbox_;
cuda::texture<float4>
  neighbour_wrapper<3>::r = neighbour_kernel::dim_<3>::r;
cuda::symbol<float3>
  neighbour_wrapper<3>::box_length = neighbour_kernel::dim_<3>::box_length;
cuda::symbol<float3>
  neighbour_wrapper<3>::cell_length = neighbour_kernel::dim_<3>::cell_length;
cuda::function<void (unsigned int*, unsigned int const*, unsigned int const*, unsigned int const*, unsigned int*)>
  neighbour_wrapper<3>::assign_cells = neighbour_kernel::assign_cells;
cuda::function<void (unsigned int*, unsigned int*)>
  neighbour_wrapper<3>::find_cell_offset = neighbour_kernel::find_cell_offset;
cuda::function<void (unsigned int*)>
  neighbour_wrapper<3>::gen_index = neighbour_kernel::gen_index;
cuda::function<void (unsigned int*, unsigned int*, unsigned int const*)>
  neighbour_wrapper<3>::update_neighbours = neighbour_kernel::update_neighbours<3>;
cuda::function<void (float4 const*, unsigned int*)>
  neighbour_wrapper<3>::compute_cell = neighbour_kernel::compute_cell<3>;

cuda::texture<float>
  neighbour_wrapper<2>::rr_cut_skin = neighbour_kernel::rr_cut_skin_;
cuda::symbol<uint2>
  neighbour_wrapper<2>::ncell = neighbour_kernel::dim_<2>::ncell;
cuda::symbol<unsigned int>
  neighbour_wrapper<2>::neighbour_size = neighbour_kernel::neighbour_size_;
cuda::symbol<unsigned int>
  neighbour_wrapper<2>::neighbour_stride = neighbour_kernel::neighbour_stride_;
cuda::symbol<unsigned int>
  neighbour_wrapper<2>::nbox = neighbour_kernel::nbox_;
cuda::texture<float4>
  neighbour_wrapper<2>::r = neighbour_kernel::dim_<2>::r;
cuda::symbol<float2>
  neighbour_wrapper<2>::box_length = neighbour_kernel::dim_<2>::box_length;
cuda::symbol<float2>
  neighbour_wrapper<2>::cell_length = neighbour_kernel::dim_<2>::cell_length;
cuda::function<void (unsigned int*, unsigned int const*, unsigned int const*, unsigned int const*, unsigned int*)>
  neighbour_wrapper<2>::assign_cells = neighbour_kernel::assign_cells;
cuda::function<void (unsigned int*, unsigned int*)>
  neighbour_wrapper<2>::find_cell_offset = neighbour_kernel::find_cell_offset;
cuda::function<void (unsigned int*)>
  neighbour_wrapper<2>::gen_index = neighbour_kernel::gen_index;
cuda::function<void (unsigned int*, unsigned int*, unsigned int const*)>
  neighbour_wrapper<2>::update_neighbours = neighbour_kernel::update_neighbours<2>;
cuda::function<void (float4 const*, unsigned int*)>
  neighbour_wrapper<2>::compute_cell = neighbour_kernel::compute_cell<2>;

}}} // namespace halmd::mdsim::gpu
