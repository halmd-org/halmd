/*
 * Copyright © 2008-2011 Peter Colberg and Felix Höfling
 * Copyright © 2012      Nicolas Höft
 * Copyright © 2016      Daniel Kirchner
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

#ifndef HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_SMOOTH_R4_KERNEL_CUH
#define HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_SMOOTH_R4_KERNEL_CUH

#include <halmd/mdsim/gpu/forces/pair_full_kernel.cuh>
#include <halmd/mdsim/gpu/forces/pair_trunc_kernel.cuh>
#include <halmd/mdsim/gpu/potentials/pair/truncations/smooth_r4_kernel.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/utility/tuple.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace pair {
namespace truncations {
namespace smooth_r4_kernel {

static texture<float4> param_;
static __constant__ float rri_smooth_;

template<typename parent_kernel>
class smooth_r4 : public parent_kernel
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
    HALMD_GPU_ENABLED smooth_r4(
        unsigned int type1, unsigned int type2
      , unsigned int ntype1, unsigned int ntype2
    )
      : parent_kernel(type1, type2, ntype1, ntype2)
      , pair_(tex1Dfetch(param_, type1 * ntype2 + type2))
    {}

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
        float_type f_abs, pot;
        tie(f_abs, pot) = parent_kernel::operator()(rr);
        pot = pot - pair_[EN_CUT];
        float_type r = sqrt(rr);
        float_type r_cut = pair_[R_CUT];
        float_type dr = r - r_cut;
        float_type x2 = dr * dr * rri_smooth_;
        float_type x4 = x2 * x2;
        float_type x4i = 1 / (1 + x4);
        // smoothing function
        float_type h0_r = x4 * x4i;
        // first derivative
        float_type h1_r = 4 * dr * rri_smooth_ * x2 * x4i * x4i;
        // apply smoothing function to obtain C¹ force function
        f_abs = h0_r * f_abs - h1_r * (pot / r);
        // apply smoothing function to obtain C² potential function
        pot = h0_r * pot;

        return make_tuple(f_abs, pot);
    }

private:
    /** adapter parameters for particle pair */
    fixed_vector<float, 3> pair_;
};

} // namespace smooth_r4_kernel

template<typename parent_kernel>
cuda::symbol<float> smooth_r4_wrapper<parent_kernel>::rri_smooth = smooth_r4_kernel::rri_smooth_;

template<typename parent_kernel>
cuda::texture<float4> smooth_r4_wrapper<parent_kernel>::param = smooth_r4_kernel::param_;

} // namespace truncations
} // namespace pair
} // namespace potentials
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_SMOOTH_R4_KERNEL_CUH */
