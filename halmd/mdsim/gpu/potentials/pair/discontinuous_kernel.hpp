/*
 * Copyright © 2016 Daniel Kirchner
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

#ifndef HALMD_MDSIM_GPU_POTENTIALS_PAIR_DISCONTINUOUS_KERNEL_HPP
#define HALMD_MDSIM_GPU_POTENTIALS_PAIR_DISCONTINUOUS_KERNEL_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace pair {
namespace discontinuous_kernel {

/**
 * indices of parameters
 */
enum {
    R_CUT       /**< cutoff length */
  , RR_CUT      /**< square of cutoff length */
  , EN_CUT      /**< potential energy at cutoff length in MD units */
};

// forward declaration for host code
template<typename parent_kernel>
class discontinuous;

} // namespace discontinuous_kernel

template<typename parent_kernel>
struct discontinuous_wrapper
{
    static cuda::texture<float4> param;
};

} // namespace pair
} // namespace potentials
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_POTENTIALS_PAIR_DISCONTINUOUS_KERNEL_HPP */
