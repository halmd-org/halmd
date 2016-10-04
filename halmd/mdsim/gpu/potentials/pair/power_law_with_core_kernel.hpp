/*
 * Copyright © 2011  Michael Kopp and Felix Höfling
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

#ifndef HALMD_MDSIM_GPU_POTENTIALS_PAIR_POWER_LAW_WITH_CORE_KERNEL_HPP
#define HALMD_MDSIM_GPU_POTENTIALS_PAIR_POWER_LAW_WITH_CORE_KERNEL_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <halmd/numeric/pow.hpp>  // std::pow is not a device function

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

template<typename float_type>
HALMD_GPU_ENABLED static inline tuple<float_type, float_type> compute(float_type const& rr
                                                                    , float_type const& sigma2
                                                                    , float_type const& epsilon
                                                                    , float_type const& core_sigma
                                                                    , unsigned short const& n)
{
    float_type rr_ss = rr / sigma2;
    // The computation of the square root can not be avoided
    // as r_core must be substracted from r but only r * r is passed.
    float_type r_s = sqrt(rr_ss);  // translates to sqrt.approx.f32 in PTX code for float_type=float (CUDA 3.2)
    float_type dri = 1 / (r_s - core_sigma);
    float_type eps_dri_n = epsilon * halmd::pow(dri, n);

    float_type en_pot = eps_dri_n;
    float_type n_eps_dri_n_1 = n * dri * eps_dri_n;
    float_type fval = n_eps_dri_n_1 / (sigma2 * r_s);

    return make_tuple(fval, en_pot);
}

} // namespace power_law_with_core_kernel

struct power_law_with_core_wrapper
{
    /** parameters for power law potential with core */
    static cuda::texture<float4> param;
};

} // namespace pair
} // namespace potentials
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_POTENTIALS_PAIR_POWER_LAW_WITH_CORE_KERNEL_HPP */
