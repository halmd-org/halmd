/*
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
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

#ifndef HALMD_MDSIM_GPU_POTENTIALS_PAIR_POWER_LAW_KERNEL_HPP
#define HALMD_MDSIM_GPU_POTENTIALS_PAIR_POWER_LAW_KERNEL_HPP

#include <halmd/numeric/blas/blas.hpp>
#include <halmd/utility/tuple.hpp>
#include <halmd/numeric/pow.hpp>  // std::pow is not a device function

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace pair {
namespace power_law_kernel {

/**
 * indices of potential parameters
 */
enum {
    EPSILON       /**< potential well depths in MD units */
  , SIGMA2        /**< square of pair separation */
  , INDEX         /**< power law index */
};

template <typename float_type>
HALMD_GPU_ENABLED static inline tuple<float_type, float_type> compute(
    float_type const& rr
  , float_type const& sigma2
  , float_type const& epsilon
  , unsigned short const& n
)
{
    float_type rri = sigma2 / rr;
    // avoid computation of square root for even powers
    float_type rni = halmd::pow(rri, n / 2);
    if (n % 2) {
        rni *= sqrt(rri); // translates to sqrt.approx.f32 in PTX code for float_type=float (CUDA 3.2)
    }
    float_type eps_rni = epsilon * rni;
    float_type fval = n * eps_rni / rr;
    float_type en_pot = eps_rni;

    return make_tuple(fval, en_pot);
}

/**
 * power law interaction potential of a pair of particles.
 *
 * @f[  U(r) = \epsilon (r/\sigma)^{-n} @f]
 */
class power_law
{
public:
    /**
     * Construct power law potential.
     */
    power_law(cudaTextureObject_t t_param) : t_param_(t_param) {}

    /**
     * Fetch potential parameters from texture cache for particle pair.
     *
     * @param type1 type of first interacting particle
     * @param type2 type of second interacting particle
     */
    HALMD_GPU_ENABLED void fetch(
        unsigned int type1, unsigned int type2
      , unsigned int ntype1, unsigned int ntype2
    );

    /**
     * Compute force and potential for interaction.
     *
     * @param rr squared distance between particles
     * @returns tuple of unit "force" @f$ -U'(r)/r @f$ and potential @f$ U(r) @f$
     *
     * @f{eqnarray*}{
     *   - U'(r) / r &=& n r^{-2} \epsilon (r/\sigma)^{-n} \\
     *   U(r) &=& \epsilon (r/\sigma)^{-n}
     * @f}
     */
    template <typename float_type>
    HALMD_GPU_ENABLED tuple<float_type, float_type> operator()(float_type rr) const
    {
        return compute(rr, pair_[SIGMA2], pair_[EPSILON], static_cast<unsigned short>(pair_[INDEX]));
    }

protected:
    /** potential parameters for particle pair */
    fixed_vector<float, 4> pair_;
    cudaTextureObject_t t_param_;
};

} // namespace power_law_kernel

struct power_law_wrapper {};

} // namespace pair
} // namespace potentials
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_POTENTIALS_PAIR_POWER_LAW_KERNEL_HPP */
