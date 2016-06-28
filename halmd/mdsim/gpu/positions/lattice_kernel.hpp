/*
 * Copyright © 2008-2010  Peter Colberg and Felix Höfling
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

#ifndef HALMD_MDSIM_GPU_POSITIONS_LATTICE_KERNEL_HPP
#define HALMD_MDSIM_GPU_POSITIONS_LATTICE_KERNEL_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace positions {

template <typename lattice_type>
struct lattice_wrapper
{
    typedef typename lattice_type::result_type vector_type;
    typedef typename lattice_type::shape_type shape_type;
    static unsigned int const dimension = vector_type::static_size;

    cuda::function<void (float4*, unsigned int, float, unsigned int, vector_type, shape_type)> lattice;

    static lattice_wrapper const kernel;
};

template <typename lattice_type>
lattice_wrapper<lattice_type> const& get_lattice_kernel()
{
    return lattice_wrapper<lattice_type>::kernel;
}

} // namespace mdsim
} // namespace gpu
} // namespace positions
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_POSITIONS_LATTICE_KERNEL_HPP */
