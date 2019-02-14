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

#ifndef HALMD_MDSIM_GPU_POTENTIALS_PAIR_MODIFIED_LENNARD_JONES_KERNEL_HPP
#define HALMD_MDSIM_GPU_POTENTIALS_PAIR_MODIFIED_LENNARD_JONES_KERNEL_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <halmd/numeric/pow.hpp>  // std::pow is not a device function

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

template<typename float_type>
HALMD_GPU_ENABLED static inline tuple<float_type, float_type> compute(float_type const& rr
                                                                    , float_type const& sigma2
                                                                    , float_type const& epsilon
                                                                    , unsigned short const& m_2
                                                                    , unsigned short const& n_2)
{
        float_type rri = sigma2 / rr;
        float_type rni = halmd::pow(rri, n_2);
        float_type rmni = (m_2 - n_2 == n_2) ? rni : halmd::pow(rri, m_2 - n_2);
        float_type eps_rni = epsilon * rni;
        float_type fval = 8 * rri * eps_rni * (m_2 * rmni - n_2) / sigma2;
        float_type en_pot = 4 * eps_rni * (rmni - 1);

        return make_tuple(fval, en_pot);

}

} // namespace modified_lennard_jones_kernel

struct modified_lennard_jones_wrapper
{
    /** Lennard-Jones potential parameters */
    static cuda::texture<float4> param;
};

} // namespace pair
} // namespace potentials
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_POTENTIALS_PAIR_MODIFIED_LENNARD_JONES_KERNEL_HPP */
