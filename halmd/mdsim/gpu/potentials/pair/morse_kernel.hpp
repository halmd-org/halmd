/*
 * Copyright © 2023       Felix Höfling
 * Copyright © 2008-2010  Peter Colberg
 * Copyright © 2020       Jaslo Ziska
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

#include <halmd/numeric/blas/blas.hpp>
#include <halmd/utility/tuple.hpp>

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
  , DISTORTION  /**< distortion factor B */
};

template <typename float_type>
HALMD_GPU_ENABLED static inline tuple<float_type, float_type> compute(
    float_type const& rr
  , float_type const& sigma
  , float_type const& epsilon
  , float_type const& r_min_sigma
  , float_type const& B     // distortion
)
{
    float_type r_sigma = sqrt(rr) / sigma / B;
    float_type exp_dr = exp(r_min_sigma / B - r_sigma);
    float_type B2 = B * B;
    float_type eps_exp_dr = epsilon * exp_dr / (2 * B2 - 1);
    float_type fval = 2 * eps_exp_dr * (exp_dr - B2) * r_sigma / rr;
    float_type en_pot = eps_exp_dr * (exp_dr - 2 * B2);

    return make_tuple(fval, en_pot);
}

/**
 * Morse potential for the interaction of a pair of particles.
 */
class morse
{
public:
    /**
     * Construct Morse's pair interaction potential.
     */
    morse(cudaTextureObject_t t_param) : t_param_(t_param) {}

    /**
     * Fetch potential parameters from texture cache for particle pair.
     *
     * @param type1 type of first interacting particle
     * @param type2 type of second interacting particle
     */
    HALMD_GPU_ENABLED void fetch_param(
        unsigned int type1, unsigned int type2
      , unsigned int ntype1, unsigned int ntype2
    );

    /**
     * Compute force and potential for interaction.
     *
     * @param rr squared distance between particles
     * @returns tuple of unit "force" @f$ -U'(r)/r @f$ and potential @f$ U(r) @f$
     */
    template <typename float_type>
    HALMD_GPU_ENABLED tuple<float_type, float_type> operator()(float_type rr) const
    {
        return morse_kernel::compute(rr, pair_[SIGMA], pair_[EPSILON], pair_[R_MIN_SIGMA], pair_[DISTORTION]);
    }

private:
    /** potential parameters for particle pair */
    fixed_vector<float, 4> pair_;
    cudaTextureObject_t t_param_;
};

} // namespace morse_kernel

struct morse_wrapper {};

} // namespace pair
} // namespace potentials
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_POTENTIALS_PAIR_MORSE_KERNEL_HPP */
