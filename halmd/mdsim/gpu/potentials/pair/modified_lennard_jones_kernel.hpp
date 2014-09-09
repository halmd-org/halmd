/*
 * Copyright © 2008-2010  Peter Colberg
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

#ifndef HALMD_MDSIM_GPU_POTENTIALS_PAIR_MODIFIED_LENNARD_JONES_KERNEL_HPP
#define HALMD_MDSIM_GPU_POTENTIALS_PAIR_MODIFIED_LENNARD_JONES_KERNEL_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace pair {
namespace modified_lennard_jones_kernel {

/**
 * indices of potential parameters
 */
enum {
    EPSILON    /**< potential well depths in MD units */
  , SIGMA2     /**< square of pair separation */
  , INDEX_M_2  /**< half-value of index of repulsion */
  , INDEX_N_2  /**< half-value of index of attraction */
};

// forward declaration for host code
class modified_lennard_jones;

} // namespace modified_lennard_jones_kernel

struct modified_lennard_jones_wrapper
{
    /** Lennard-Jones potential parameters */
    static cuda::texture<float4> param;
    /** squared cutoff radius and energy shift */
    static cuda::texture<float2> rr_en_cut;
};

} // namespace pair
} // namespace potentials
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_POTENTIALS_PAIR_MODIFIED_LENNARD_JONES_KERNEL_HPP */
