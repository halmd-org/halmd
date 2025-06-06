/*
 * Copyright © 2025 Giorgia Marcelli
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

#ifndef HALMD_MDSIM_GPU_VELOCITIES_RESCALE_KERNEL_HPP
#define HALMD_MDSIM_GPU_VELOCITIES_RESCALE_KERNEL_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <halmd/mdsim/type_traits.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace velocities {

// return codes of CUDA kernel
enum {
    success = 0x0   // everything OK
  , warning = 0x1   // zero velocity encountered, workaround used
  , failure = 0x2   // target energy too low, terminate
};

template <int dimension, typename float_type>
struct rescale_wrapper
{
    typedef typename type_traits<dimension, float>::gpu::coalesced_vector_type coalesced_vector_type;
    typedef typename type_traits<4, float_type>::gpu::ptr_type ptr_type;

    // declare the rescale kernel function
    cuda::function<void(ptr_type, float const*, uint, float, int*)> rescale;

    static rescale_wrapper kernel;
};

} // namespace velocities
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif // HALMD_MDSIM_GPU_VELOCITIES_RESCALE_KERNEL_HPP
