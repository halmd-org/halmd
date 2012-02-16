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

#ifndef HALMD_MDSIM_GPU_POTENTIALS_MORSE_KERNEL_HPP
#define HALMD_MDSIM_GPU_POTENTIALS_MORSE_KERNEL_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace morse_kernel {

/**
 * indices of potential parameters in float4 array
 */
enum {
    EPSILON     /**< depth of potential well in MD units */
  , SIGMA       /**< width of potential well in MD units */
  , R_MIN_SIGMA /**< position of potential well in units of sigma */
  , EN_CUT      /**< potential energy at cutoff radius in MD units */
};

// forward declaration for host code
class morse;

} // namespace morse_kernel

struct morse_wrapper
{
    /** potential parameters */
    static cuda::texture<float4> param;
    /** squared cutoff radius */
    static cuda::texture<float> rr_cut;
};

} // namespace potentials
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_POTENTIALS_MORSE_KERNEL_HPP */
