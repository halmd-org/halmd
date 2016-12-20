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

#ifndef HALMD_MDSIM_GPU_POTENTIALS_PAIR_MORSE_KERNEL_HPP
#define HALMD_MDSIM_GPU_POTENTIALS_PAIR_MORSE_KERNEL_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace pair {
namespace morse_kernel {

/**
 * indices of potential parameters in float4 array
 */
enum {
    EPSILON     /**< depth of potential well in MD units */
  , SIGMA       /**< width of potential well in MD units */
  , R_MIN_SIGMA /**< position of potential well in units of sigma */
};

// forward declaration for host code
class morse;

template<typename float_type>
HALMD_GPU_ENABLED static inline tuple<float_type, float_type> compute(float_type const& rr
                                                                    , float_type const& sigma
                                                                    , float_type const& epsilon
                                                                    , float_type const& r_min)
{
    float_type r_sigma = sqrt(rr) / sigma;
    float_type exp_dr = exp(r_min - r_sigma);
    float_type eps_exp_dr = epsilon * exp_dr;
    float_type fval = 2 * eps_exp_dr * (exp_dr - 1) * r_sigma / rr;
    float_type en_pot = eps_exp_dr * (exp_dr - 2);

    return make_tuple(fval, en_pot);
}

} // namespace morse_kernel

struct morse_wrapper
{
    /** potential parameters */
    static cuda::texture<float4> param;
};

} // namespace pair
} // namespace potentials
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_POTENTIALS_PAIR_MORSE_KERNEL_HPP */
