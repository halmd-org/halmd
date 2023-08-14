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

#ifndef HALMD_MDSIM_GPU_POTENTIALS_PAIR_LENNARD_JONES_SIMPLE_KERNEL_HPP
#define HALMD_MDSIM_GPU_POTENTIALS_PAIR_LENNARD_JONES_SIMPLE_KERNEL_HPP

#include <halmd/utility/tuple.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace pair {
namespace lennard_jones_simple_kernel {

template <typename float_type>
HALMD_GPU_ENABLED static inline tuple<float_type, float_type> compute(float_type rr)
{
    const float_type sigma2 = 1;
    const float_type epsilon = 1;
    float_type rri = sigma2 / rr;
    float_type ri6 = rri * rri * rri;
    float_type eps_ri6 = epsilon * ri6;
    float_type fval = 48 * rri * eps_ri6 * (ri6 - 0.5f) / sigma2;
    float_type en_pot = 4 * eps_ri6 * (ri6 - 1);

    return make_tuple(fval, en_pot);
}

/**
 * Lennard-Jones interaction for a simple fluid of a single species.
 */
class lennard_jones_simple
{
public:
    /**
     * @param type1 type of first interacting particle
     * @param type2 type of second interacting particle
     */
    HALMD_GPU_ENABLED void fetch_param(
        unsigned int type1, unsigned int type2
      , unsigned int ntype1, unsigned int ntype2
    ) {}

    /**
     * Compute force and potential for interaction.
     *
     * @param rr squared distance between particles
     * @returns tuple of unit "force" @f$ -U'(r)/r @f$ and potential @f$ U(r) @f$
     */
    template <typename float_type>
    HALMD_GPU_ENABLED tuple<float_type, float_type> operator()(float_type rr) const
    {
        return lennard_jones_simple_kernel::compute(rr);
    }
};

} // namespace lennard_jones_kernel
} // namespace pair
} // namespace potentials
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_POTENTIALS_PAIR_LENNARD_JONES_KERNEL_HPP */
