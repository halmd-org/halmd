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

#ifndef HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_SHIFTED_KERNEL_CUH
#define HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_SHIFTED_KERNEL_CUH

#include <halmd/mdsim/gpu/forces/pair_full_kernel.cuh>
#include <halmd/mdsim/gpu/forces/pair_trunc_kernel.cuh>
#include <halmd/mdsim/gpu/potentials/pair/truncations/shifted_kernel.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/utility/tuple.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace pair {
namespace truncations {
namespace shifted_kernel {

static texture<float2> param_;

template<typename parent_kernel>
class shifted
  : public parent_kernel
{
public:
    /**
     * Construct Smoothing Function.
     *
     * Fetch potential parameters from texture cache for particle pair.
     *
     * @param type1 type of first interacting particle
     * @param type2 type of second interacting particle
     */
    HALMD_GPU_ENABLED shifted(
        unsigned int type1, unsigned int type2
      , unsigned int ntype1, unsigned int ntype2
    )
      : parent_kernel(type1, type2, ntype1, ntype2)
      , pair_(tex1Dfetch(param_, type1 * ntype2 + type2))
    {
    }

    /**
     * Check whether particles are in interaction range.
     *
     * @param rr squared distance between particles
     */
    template <typename float_type>
    HALMD_GPU_ENABLED bool within_range(float_type rr) const
    {
        return (rr < pair_[RR_CUT]);
    }

    /**
     * Compute force and potential for interaction.
     *
     * @param rr squared distance between particles
     * @returns tuple of unit "force" @f$ -U'(r)/r @f$ and potential @f$ U(r) @f$
     */
    template <typename float_type>
    HALMD_GPU_ENABLED tuple<float_type, float_type> operator()(float_type rr) const
    {
        float_type f_abs, en_pot;
        tie(f_abs, en_pot) = parent_kernel::operator()(rr);
        en_pot = en_pot - pair_[EN_CUT];
        return make_tuple(f_abs, en_pot);
    }

private:
    /** adapter parameters for particle pair */
    fixed_vector<float, 2> pair_;
};

} // namespace shifted_kernel

template<typename parent_kernel>
cuda::texture<float2> shifted_wrapper<parent_kernel>::param = shifted_kernel::param_;

} // namespace truncations
} // namespace pair
} // namespace potentials
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_SHIFTED_KERNEL_CUH */
