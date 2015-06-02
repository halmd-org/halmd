/*
 * Copyright © 2014-2015 Sutapa Roy
 * Copyright © 2014-2015 Felix Höfling
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

#ifndef HALMD_MDSIM_GPU_POTENTIALS_EXTERNAL_SLIT_KERNEL_HPP
#define HALMD_MDSIM_GPU_POTENTIALS_EXTERNAL_SLIT_KERNEL_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace external {
namespace slit_kernel {

/**
 * indices of potential parameters
 */
enum {
    EPSILON     /* interaction strength */
  , SIGMA       /* interaction range */
  , WETTING     /* wetting parameter */
  , CUTOFF      /* cutoff distance */
};

// forward declaration for host code
template <int dimension>
class slit;

} // namespace slit_kernel

struct slit_wrapper
{
    /** potential parameters */
    static cuda::texture<float4> param_geometry;
    static cuda::texture<float4> param_potential;
    static cuda::symbol<float> smoothing;
    static cuda::symbol<int> nwall;
};

} // namespace external
} // namespace potentials
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_POTENTIALS_EXTERNAL_SLIT_KERNEL_HPP */
