/*
 * Copyright © 2016       Manuel Dibak
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

#ifndef HALMD_MDSIM_GPU_ORIENTATIONS_UNIFORM_KERNEL_HPP
#define HALMD_MDSIM_GPU_ORIENTATIONS_UNIFORM_KERNEL_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <halmd/mdsim/type_traits.hpp>
#include <halmd/numeric/blas/fixed_vector.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace orientations {

template <int dimension, typename float_type, typename rng_type>
struct uniform_wrapper
{
    typedef typename type_traits<4, float_type>::gpu::ptr_type ptr_type;

    cuda::function<void (ptr_type, unsigned int, rng_type)> uniform;

    static uniform_wrapper const kernel;
};

} // namespace mdsim
} // namespace gpu
} // namespace orientations
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_ORIENTATIONS_UNIFORM_KERNEL_HPP */
