/* Lennard-Jones fluid simulation using CUDA
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

#ifndef LJGPU_MDSIM_IMPL_HPP
#define LJGPU_MDSIM_IMPL_HPP

namespace ljgpu
{

template <int dimension>
class mdsim_impl;

template <int dimension>
class ljfluid_impl_base;

template <int dimension>
class ljfluid_impl_gpu_base;

template <int dimension>
class ljfluid_impl_gpu_square;

template <int dimension>
class ljfluid_impl_gpu_cell;

template <int dimension>
class ljfluid_impl_gpu_neighbour;

template <int dimension>
class ljfluid_impl_host;

template <int dimension>
class hardsphere_impl;

} // namespace ljgpu

#endif /* ! LJGPU_MDSIM_IMPL_HPP */
