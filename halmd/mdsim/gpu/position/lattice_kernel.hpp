/*
 * Copyright © 2008-2010  Peter Colberg and Felix Höfling
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

#ifndef HALMD_MDSIM_GPU_POSITION_LATTICE_KERNEL_HPP
#define HALMD_MDSIM_GPU_POSITION_LATTICE_KERNEL_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>

namespace halmd
{
namespace mdsim { namespace gpu { namespace position
{

template <int dimension>
struct lattice_wrapper
{
    cuda::function<void (float4*, uint, float)> fcc;
    cuda::function<void (float4*, uint, float)> sc;
    static lattice_wrapper const kernel;
};

template <int dimension>
lattice_wrapper<dimension> const& get_lattice_kernel()
{
    return lattice_wrapper<dimension>::kernel;
}

}}} // namespace mdsim::gpu::position

} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_POSITION_LATTICE_KERNEL_HPP */
