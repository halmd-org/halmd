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

#ifndef HALMD_MDSIM_GPU_POSITIONS_LATTICE_KERNEL_HPP
#define HALMD_MDSIM_GPU_POSITIONS_LATTICE_KERNEL_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <halmd/mdsim/type_traits.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace positions {

template <int dimension>
struct lattice_wrapper
{
    typedef typename type_traits<dimension, float>::gpu::vector_type vector_type;
    typedef typename type_traits<dimension, unsigned int>::gpu::vector_type index_type;

    /** offset: origin of particle lattice */
    cuda::symbol<vector_type> offset;
    /** number of cells per dimension */
    cuda::symbol<index_type> ncell;

    cuda::function<void (float4*, uint, float, uint)> fcc;
    cuda::function<void (float4*, uint, float, uint)> sc;
    static lattice_wrapper const kernel;
};

template <int dimension>
lattice_wrapper<dimension> const& get_lattice_kernel()
{
    return lattice_wrapper<dimension>::kernel;
}

} // namespace mdsim
} // namespace gpu
} // namespace positions
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_POSITIONS_LATTICE_KERNEL_HPP */
