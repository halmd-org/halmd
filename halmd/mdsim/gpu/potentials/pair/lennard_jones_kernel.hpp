/*
 * Copyright © 2008-2010  Peter Colberg
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

#ifndef HALMD_MDSIM_GPU_POTENTIALS_PAIR_LENNARD_JONES_KERNEL_HPP
#define HALMD_MDSIM_GPU_POTENTIALS_PAIR_LENNARD_JONES_KERNEL_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace pair {
namespace lennard_jones_kernel {

/**
 * indices of potential parameters
 */
enum {
    RR_CUT      /**< square of cutoff length */
  , EPSILON     /**< potential well depths in MD units */
  , SIGMA2      /**< square of pair separation */
  , EN_CUT      /**< potential energy at cutoff length in MD units */
};

// forward declaration for host code
class lennard_jones;

} // namespace lennard_jones_kernel

struct lennard_jones_wrapper
{
    /** Lennard-Jones potential parameters */
    static cuda::texture<float4> param;
};

} // namespace pair
} // namespace potentials
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_POTENTIALS_PAIR_LENNARD_JONES_KERNEL_HPP */
