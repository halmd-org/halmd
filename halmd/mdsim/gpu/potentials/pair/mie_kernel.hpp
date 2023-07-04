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

#ifndef HALMD_MDSIM_GPU_POTENTIALS_PAIR_MODIFIED_LENNARD_JONES_KERNEL_HPP
#define HALMD_MDSIM_GPU_POTENTIALS_PAIR_MODIFIED_LENNARD_JONES_KERNEL_HPP

#include <halmd/numeric/blas/blas.hpp>
#include <halmd/numeric/pow.hpp>  // std::pow is not a device function
#include <halmd/utility/tuple.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace pair {
namespace mie_kernel {

/**
 * indices of potential parameters
 */
enum {
    EPSILON_C  /**< potential well depths in MD units */
  , SIGMA2     /**< square of pair separation */
  , INDEX_M_2  /**< half-value of index of repulsion */
  , INDEX_N_2  /**< half-value of index of attraction */
};

template<typename float_type>
HALMD_GPU_ENABLED static inline tuple<float_type, float_type> compute(
    float_type const& rr
  , float_type const& sigma2
  , float_type const& epsilon_C
  , unsigned short const& m_2
  , unsigned short const& n_2)
{
    float_type rri = sigma2 / rr;
    float_type rni = halmd::pow(rri, n_2);
    float_type rmni = (m_2 - n_2 == n_2) ? rni : halmd::pow(rri, m_2 - n_2);
    float_type eps_rni = epsilon_C * rni;
    float_type fval = 2 * rri * eps_rni * (m_2 * rmni - n_2) / sigma2;
    float_type en_pot = eps_rni * (rmni - 1);

    return make_tuple(fval, en_pot);
}

/**
 *
 */
class mie
{
public:
    /**
     * Construct Lennard-Jones pair interaction potential.
     */
    mie(cudaTextureObject_t t_param) : t_param_(t_param) {}

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
     *
     * @f{eqnarray*}{
     *   - U'(r) / r &=& C(m, n) r^{-2} \epsilon (\sigma/r)^{n} \left[ m (\sigma/r)^{m-n} - n \right] \\
     *   U(r) &=& C(m, n) \epsilon (\sigma/r)^{n} \left[ (\sigma/r)^{m-n} - 1 \right]
     * @f}
     */
    template <typename float_type>
    HALMD_GPU_ENABLED tuple<float_type, float_type> operator()(float_type rr) const
    {
        return compute(rr, pair_[SIGMA2], pair_[EPSILON_C]
          , static_cast<unsigned short>(pair_[INDEX_M_2])
          , static_cast<unsigned short>(pair_[INDEX_N_2])
        );
    }

private:
    /** potential parameters for particle pair */
    fixed_vector<float, 4> pair_;
    cudaTextureObject_t t_param_;
};

} // namespace mie_kernel

struct mie_wrapper {};

} // namespace pair
} // namespace potentials
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_POTENTIALS_PAIR_MODIFIED_LENNARD_JONES_KERNEL_HPP */
