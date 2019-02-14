/*
 * Copyright © 2016 Daniel Kirchner
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

#ifndef HALMD_MDSIM_GPU_POTENTIALS_PAIR_POWER_LAW_HARD_CORE_KERNEL_CUH
#define HALMD_MDSIM_GPU_POTENTIALS_PAIR_POWER_LAW_HARD_CORE_KERNEL_CUH

#include <halmd/mdsim/gpu/potentials/pair/adapters/hard_core_kernel.cuh>
#include <halmd/mdsim/gpu/potentials/pair/power_law_kernel.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace pair {
namespace adapters {

namespace hard_core_kernel {

/** explicit instantiation for power_law kernel */
template <>
class hard_core<power_law_kernel::power_law> : public power_law_kernel::power_law
{
public:
    /**
     * Construct Hard Core Adapter.
     *
     * Fetch core parameter from texture cache for particle pair.
     *
     * @param type1 type of first interacting particle
     * @param type2 type of second interacting particle
     */
    HALMD_GPU_ENABLED hard_core(
        unsigned int type1, unsigned int type2
      , unsigned int ntype1, unsigned int ntype2
    )
      : power_law_kernel::power_law(type1, type2, ntype1, ntype2)
      , core_sigma_(tex1Dfetch(param_, type1 * ntype2 + type2))
    {}

    /**
     * Compute force and potential for interaction.
     *
     * @param rr squared distance between particles
     * @returns tuple of unit "force" @f$ -U'(r)/r @f$ and potential @f$ U(r) @f$
     */
    template <typename float_type>
    HALMD_GPU_ENABLED tuple<float_type, float_type> operator()(float_type rr) const
    {
        unsigned short n = static_cast<unsigned short>(pair_[power_law_kernel::INDEX]);
        float_type rr_ss = rr / pair_[power_law_kernel::SIGMA2];
        // The computation of the square root can not be avoided
        // as r_core must be substracted from r but only r * r is passed.
        float_type r_s = sqrtf(rr_ss);  // translates to sqrt.approx.f32 in PTX code for float_type=float (CUDA 3.2)
        float_type dri = 1 / (r_s - core_sigma_);
        float_type eps_dri_n = pair_[power_law_kernel::EPSILON] * halmd::pow(dri, n);

        float_type en_pot = eps_dri_n;
        float_type n_eps_dri_n_1 = n * dri * eps_dri_n;
        float_type fval = n_eps_dri_n_1 / (pair_[power_law_kernel::SIGMA2] * r_s);

        return make_tuple(fval, en_pot);
    }

private:
    float core_sigma_;
};

} // namespace hard_core_kernel
} // namespace adapters
} // namespace pair
} // namespace potentials
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_POTENTIALS_PAIR_POWER_LAW_HARD_CORE_KERNEL_CUH */
