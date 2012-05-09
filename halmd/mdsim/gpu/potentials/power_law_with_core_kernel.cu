/*
 * Copyright © 2011-2012  Michael Kopp and Felix Höfling
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
#include <halmd/mdsim/gpu/potentials/power_law_with_core_kernel.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/numeric/pow.hpp>  // std::pow is not a device function
#include <halmd/utility/tuple.hpp>
#include <halmd/mdsim/smoothers/localr4.hpp>
#include <halmd/mdsim/smoothers/nosmooth.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace power_law_with_core_kernel {

/** array of potential parameters for all combinations of particle types */
static texture<float4> param_;
/** squares of potential cutoff radius and energy shift for all combinations of particle types */
static texture<float2> rr_en_cut_;

/**
 * power law interaction potential of a pair of particles.
 *
 * @f[  U(r) = \epsilon \left(\frac{\sigma}{r - r_\mathrm{core}}\right)^n @f]
 */
class power_law_with_core
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
    HALMD_GPU_ENABLED power_law_with_core(
        unsigned int type1, unsigned int type2
      , unsigned int ntype1, unsigned int ntype2
    )
      : pair_(tex1Dfetch(param_, type1 * ntype2 + type2))
      , pair_rr_en_cut_(tex1Dfetch(rr_en_cut_, type1 * ntype2 + type2))
    {}

    /**
     * Returns square of cutoff distance.
     */
    HALMD_GPU_ENABLED float rr_cut() const
    {
        return pair_rr_en_cut_[0];
    }

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
     *   U(r) &=& \epsilon \left(\frac{r - r_\mathrm{core}}{\sigma}\right)^{-n} \\
     *   - \frac{U'(r)}{r} &=& n \frac{n}{r(r-r_\mathrm{core})} U(r) \\
     *   r \partial_r r \partial_r U(r) &=& n \frac{r (n r + r_\mathrm{core}) }{ (r-r_\mathrm{core})^2 } U(r)
     * @f}
     */
    template <typename float_type>
    HALMD_GPU_ENABLED tuple<float_type, float_type, float_type> operator()(float_type rr) const
    {
        float_type rr_ss = rr / pair_[SIGMA2];
        // The computation of the square root can not be avoided
        // as r_core must be substracted from r but only r * r is passed.
        float_type r_s = sqrt(rr_ss);  // translates to sqrt.approx.f32 in PTX code for float_type=float (CUDA 3.2)
        float_type dri = 1 / (r_s - pair_[CORE_SIGMA]);
        unsigned short n = static_cast<unsigned short>(pair_[INDEX]);
        float_type eps_dri_n = pair_[EPSILON] * halmd::pow(dri, n);

        float_type en_pot = eps_dri_n - pair_rr_en_cut_[1];
        float_type n_eps_dri_n_1 = n * dri * eps_dri_n;
        float_type fval = n_eps_dri_n_1 / (pair_[SIGMA2] * r_s);
        float_type hvir = n_eps_dri_n_1 * ((n + 1) * dri * rr_ss - r_s);

        return make_tuple(fval, en_pot, hvir);
    }

private:
    /** potential parameters for particle pair */
    fixed_vector<float, 4> pair_;
    /** squared cutoff radius and energy shift for particle pair */
    fixed_vector<float, 2> pair_rr_en_cut_;
};

} // namespace power_law_with_core_kernel

cuda::texture<float4> power_law_with_core_wrapper::param = power_law_with_core_kernel::param_;
cuda::texture<float2> power_law_with_core_wrapper::rr_en_cut = power_law_with_core_kernel::rr_en_cut_;

} // namespace potentials

// explicit instantiation of force kernels
namespace forces {

using namespace halmd::mdsim::gpu::potentials::power_law_with_core_kernel;
using namespace halmd::mdsim::smoothers;

template class pair_full_wrapper<3, power_law_with_core>;
template class pair_full_wrapper<2, power_law_with_core>;

template class pair_trunc_wrapper<3, power_law_with_core, nosmooth>;
template class pair_trunc_wrapper<2, power_law_with_core, nosmooth>;
template class pair_trunc_wrapper<3, power_law_with_core, localr4<float> >;
template class pair_trunc_wrapper<2, power_law_with_core, localr4<float> >;

} // namespace forces

} // namespace gpu
} // namespace mdsim
} // namespace halmd
