/* Lennard-Jones fluid kernel
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

#ifndef LJGPU_LJFLUID_GPU_LATTICE_HPP
#define LJGPU_LJFLUID_GPU_LATTICE_HPP

#include <cuda_wrapper.hpp>

namespace ljgpu { namespace gpu { namespace lattice
{

extern cuda::function<void (float2*, uint, float), void (float4*, uint, float)> fcc;
extern cuda::function<void (float2*, uint, float), void (float4*, uint, float)> sc;

}}} // namespace ljgpu::gpu::lattice

#endif /* ! LJGPU_LJFLUID_GPU_LATTICE_HPP */
