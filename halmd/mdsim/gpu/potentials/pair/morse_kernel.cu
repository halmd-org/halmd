/*
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
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

#include <halmd/mdsim/gpu/forces/pair_full_kernel.cuh>
#include <halmd/mdsim/gpu/forces/pair_trunc_kernel.cuh>
#include <halmd/mdsim/gpu/potentials/pair/local_r4_kernel.cuh>
#include <halmd/mdsim/gpu/potentials/pair/morse_kernel.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/utility/tuple.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace pair {
namespace morse_kernel {

/** array of potential parameters for all combinations of particle types */
static texture<float4> param_;
/** squares of potential cutoff radius for all combinations of particle types */
static texture<float> rr_cut_;

/**
 * Morse potential for the interaction of a pair of particles.
 */
class morse
{
public:
    /**
     * Construct Morse's pair interaction potential.
     *
     * Fetch potential parameters from texture cache for particle pair.
     *
     * @param type1 type of first interacting particle
     * @param type2 type of second interacting particle
     */
    HALMD_GPU_ENABLED morse(
        unsigned int type1, unsigned int type2
      , unsigned int ntype1, unsigned int ntype2
    )
      : pair_(tex1Dfetch(param_, type1 * ntype2 + type2))
      , pair_rr_cut_(tex1Dfetch(rr_cut_, type1 * ntype2 + type2))
    {}

    /**
     * Returns square of cutoff distance.
     */
    HALMD_GPU_ENABLED float rr_cut() const
    {
        return pair_rr_cut_;
    }

    /**
     * Check whether particles are in interaction range.
     *
     * @param rr squared distance between particles
     */
    template <typename float_type>
    HALMD_GPU_ENABLED bool within_range(float_type rr) const
    {
        return (rr < pair_rr_cut_);
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
        float_type r_sigma = sqrt(rr) / pair_[SIGMA];
        float_type exp_dr = exp(pair_[R_MIN_SIGMA] - r_sigma);
        float_type eps_exp_dr = pair_[EPSILON] * exp_dr;
        float_type fval = 2 * eps_exp_dr * (exp_dr - 1) * r_sigma / rr;
        float_type en_pot = eps_exp_dr * (exp_dr - 2) - pair_[EN_CUT];

        return make_tuple(fval, en_pot);
    }

private:
    /** potential parameters for particle pair */
    fixed_vector<float, 4> pair_;
    /** squared cutoff radius for particle pair */
    float pair_rr_cut_;
};

} // namespace morse_kernel

cuda::texture<float4> morse_wrapper::param = morse_kernel::param_;
cuda::texture<float> morse_wrapper::rr_cut = morse_kernel::rr_cut_;
template class local_r4_wrapper<morse_kernel::morse>;

} // namespace pair
} // namespace potentials

// explicit instantiation of force kernels
namespace forces {

using namespace halmd::mdsim::gpu::potentials::pair::morse_kernel;
using namespace halmd::mdsim::gpu::potentials::pair::local_r4_kernel;

template class pair_full_wrapper<3, morse>;
template class pair_full_wrapper<2, morse>;

template class pair_trunc_wrapper<3, morse>;
template class pair_trunc_wrapper<2, morse>;
template class pair_trunc_wrapper<3, local_r4<morse> >;
template class pair_trunc_wrapper<2, local_r4<morse> >;

} // namespace forces

} // namespace gpu
} // namespace mdsim
} // namespace halmd
