/*
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

#ifndef HALMD_MDSIM_GPU_POTENTIALS_PAIR_LENNARD_JONES_KERNEL_HPP
#define HALMD_MDSIM_GPU_POTENTIALS_PAIR_LENNARD_JONES_KERNEL_HPP

#include <halmd/utility/tuple.hpp>
#include <halmd/numeric/blas/blas.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace pair {
namespace lennard_jones_kernel {

/**
 * indices of potential parameters
 */
enum {
    EPSILON     /**< potential well depths in MD units */
  , SIGMA2      /**< square of pair separation */
};

template <typename float_type>
HALMD_GPU_ENABLED static inline tuple<float_type, float_type> compute(
    float_type const& rr
  , float_type const& sigma2
  , float_type const& epsilon
)
{
    float_type rri =  sigma2 / rr;
    float_type ri6 = rri * rri * rri;
    float_type eps_ri6 = epsilon * ri6;
    float_type fval = 48 * rri * eps_ri6 * (ri6 - 0.5f) / sigma2;
    float_type en_pot = 4 * eps_ri6 * (ri6 - 1);
    return make_tuple(fval, en_pot);
}

/**
 * Lennard-Jones interaction of a pair of particles.
 */
class lennard_jones
{
public:
    /**
     * Construct Lennard-Jones pair interaction potential.
     */
    lennard_jones(cudaTextureObject_t t_param) : t_param_(t_param) {}

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
     */
    template <typename float_type>
    HALMD_GPU_ENABLED tuple<float_type, float_type> operator()(float_type rr) const;

private:
    /** potential parameters for particle pair */
    fixed_vector<float, 2> pair_;
    cudaTextureObject_t t_param_;
};

template <typename float_type>
HALMD_GPU_ENABLED tuple<float_type, float_type> lennard_jones::operator()(float_type rr) const
{
    return compute(rr, pair_[SIGMA2], pair_[EPSILON]);
}

} // namespace lennard_jones_kernel

struct lennard_jones_wrapper {};

} // namespace pair
} // namespace potentials
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_POTENTIALS_PAIR_LENNARD_JONES_KERNEL_HPP */
