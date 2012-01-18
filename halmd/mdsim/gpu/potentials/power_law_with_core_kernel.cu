/*
 * Copyright © 2011  Michael Kopp
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

#include <halmd/algorithm/gpu/tuple.cuh>
#include <halmd/mdsim/gpu/forces/pair_full_kernel.cuh>
#include <halmd/mdsim/gpu/forces/pair_trunc_kernel.cuh>
#include <halmd/mdsim/gpu/potentials/power_law_with_core_kernel.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/numeric/pow.hpp>  // std::pow is not a device function
#include <halmd/utility/gpu/variant.cuh>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace power_law_with_core_kernel {

using algorithm::gpu::tuple;
using algorithm::gpu::make_tuple;

/** array of potential parameters for all combinations of particle types */
static texture<float4> param_;
/** squares of potential cutoff radius and energy shift for all combinations of particle types */
static texture<float2> rr_en_cut_;

/**
 * power law interaction potential of a pair of particles.
 *
 * @f[  U(r) = \epsilon \left(\frac{\sigma}{r - r_\text{core}}\right)^n @f]
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
    HALMD_GPU_ENABLED power_law_with_core(unsigned int type1, unsigned int type2)
      : pair_(
            tex1Dfetch(param_, symmetric_matrix::lower_index(type1, type2))
        )
      , pair_rr_en_cut_(
            tex1Dfetch(rr_en_cut_, symmetric_matrix::lower_index(type1, type2))
        ) {}

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
     *   - \frac{U'(r)}{r} &=& n \epsilon \frac{\sigma^n}{r (r-r_\text{core})^{n+1}} = \frac{n}{r(r-r_\text{core})} U(r) \\
     *   U(r) &=& \epsilon \left(\frac{\sigma}{r - r_\text{core}}\right)^n \\
     *   r \partial_r r \partial_r U(r) &=& n U(r) \left( \frac{r (r_\text{core} + n r) }{ (r-r_\text{core})^2 } \right)
     * @f}
     */
    template <typename float_type>
    HALMD_GPU_ENABLED tuple<float_type, float_type, float_type> operator()(float_type rr) const
    {
        float_type rrs = rr / pair_[SIGMA2];
        // It's not possible to avoid computation of a squareroot as r_core
        // must be substracted from r but only r*r is passed.
        float_type rs = sqrt(rrs);
        unsigned short n = static_cast<unsigned short>(pair_[INDEX]);
        float coreS = pair_[CORE_SIGMA];
        float_type en_pot = pair_[EPSILON] * halmd::pow(rs - coreS, -n);

        float_type fval = en_pot * n / (pair_[SIGMA2] * rs * (rs - coreS));
        float_type hvir = fval * rr * (coreS + n*rs)/(rs - coreS);

        en_pot -= pair_rr_en_cut_[1];

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

using potentials::power_law_with_core_kernel::power_law_with_core;

template class pair_full_wrapper<3, power_law_with_core>;
template class pair_full_wrapper<2, power_law_with_core>;

template class pair_trunc_wrapper<3, power_law_with_core>;
template class pair_trunc_wrapper<2, power_law_with_core>;

} // namespace forces

} // namespace gpu
} // namespace mdsim
} // namespace halmd
