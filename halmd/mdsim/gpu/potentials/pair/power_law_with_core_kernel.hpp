/*
 * Copyright © 2011  Michael Kopp and Felix Höfling
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

#ifndef HALMD_MDSIM_GPU_POTENTIALS_PAIR_POWER_LAW_WITH_CORE_KERNEL_HPP
#define HALMD_MDSIM_GPU_POTENTIALS_PAIR_POWER_LAW_WITH_CORE_KERNEL_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace pair {
namespace power_law_with_core_kernel {

/**
 * indices of potential parameters
 */
enum {
    EPSILON       /**< potential well depths in MD units */
  , SIGMA2        /**< square of pair separation */
  , INDEX         /**< power law index */
  , CORE_SIGMA    /**< core radius in units of sigma */
};

// forward declaration for host code
class power_law_with_core;

} // namespace power_law_with_core_kernel

struct power_law_with_core_wrapper
{
    /** parameters for power law potential with core */
    static cuda::texture<float4> param;
    /** squared cutoff radius and energy shift */
    static cuda::texture<float2> rr_en_cut;
};

} // namespace pair
} // namespace potentials
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_POTENTIALS_PAIR_POWER_LAW_WITH_CORE_KERNEL_HPP */
