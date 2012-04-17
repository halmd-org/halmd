/*
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <halmd/mdsim/gpu/forces/pair_full_kernel.cuh>
#include <halmd/mdsim/gpu/forces/pair_trunc_kernel.cuh>
#include <halmd/mdsim/gpu/potentials/power_law_kernel.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/numeric/pow.hpp>  // std::pow is not a device function
#include <halmd/utility/tuple.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace power_law_kernel {

/** array of potential parameters for all combinations of particle types */
static texture<float4> param_;
/** squares of potential cutoff radius and energy shift for all combinations of particle types */
static texture<float2> rr_en_cut_;

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
     *
     * Fetch potential parameters from texture cache for particle pair.
     *
     * @param type1 type of first interacting particle
     * @param type2 type of second interacting particle
     */
    HALMD_GPU_ENABLED power_law(
        unsigned int type1, unsigned int type2
      , unsigned int ntype1, unsigned int ntype2
    )
      : pair_(tex1Dfetch(param_, type1 * ntype2 + type2))
      , pair_rr_en_cut_(tex1Dfetch(rr_en_cut_, type1 * ntype2 + type2))
    {}

    /**
     * Check whether particles are in interaction range.
     *
     * @param rr squared distance between particles
     */
    template <typename float_type>
    HALMD_GPU_ENABLED bool within_range(float_type rr) const
    {
        return (rr < pair_rr_en_cut_[0]);
    }

    /**
     * Compute force and potential for interaction.
     *
     * @param rr squared distance between particles
     * @returns tuple of unit "force" @f$ -U'(r)/r @f$, potential @f$ U(r) @f$,
     * and hypervirial @f$ r \partial_r r \partial_r U(r) @f$
     *
     * @f{eqnarray*}{
     *   - U'(r) / r &=& n r^{-2} \epsilon (r/\sigma)^{-n} \\
     *   U(r) &=& \epsilon (r/\sigma)^{-n} \\
     *   r \partial_r r \partial_r U(r) &=& n^2 \epsilon (r/\sigma)^{-n}
     * @f}
     */
    template <typename float_type>
    HALMD_GPU_ENABLED tuple<float_type, float_type, float_type> operator()(float_type rr) const
    {
        float_type rri = pair_[SIGMA2] / rr;
        unsigned short n = static_cast<unsigned short>(pair_[INDEX]);
        // avoid computation of square root for even powers
        float_type rni = halmd::pow(rri, n / 2);
        if (n % 2) {
            rni *= sqrt(rri); // translates to sqrt.approx.f32 in PTX code for float_type=float (CUDA 3.2)
        }
        float_type eps_rni = pair_[EPSILON] * rni;
        float_type fval = n * eps_rni / rr;
        float_type en_pot = eps_rni - pair_rr_en_cut_[1];
        float_type hvir = n * n * eps_rni;

        return make_tuple(fval, en_pot, hvir);
    }

private:
    /** potential parameters for particle pair */
    fixed_vector<float, 4> pair_;
    /** squared cutoff radius and energy shift for particle pair */
    fixed_vector<float, 2> pair_rr_en_cut_;
};

} // namespace power_law_kernel

cuda::texture<float4> power_law_wrapper::param = power_law_kernel::param_;
cuda::texture<float2> power_law_wrapper::rr_en_cut = power_law_kernel::rr_en_cut_;

} // namespace potentials

// explicit instantiation of force kernels
namespace forces {

using potentials::power_law_kernel::power_law;

template class pair_full_wrapper<3, power_law>;
template class pair_full_wrapper<2, power_law>;

template class pair_trunc_wrapper<3, power_law>;
template class pair_trunc_wrapper<2, power_law>;

} // namespace forces

} // namespace gpu
} // namespace mdsim
} // namespace halmd
