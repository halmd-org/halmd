/*
 * Copyright © 2008-2010  Peter Colberg and Felix Höfling
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

#include <halmd/mdsim/gpu/forces/lj_kernel.hpp>
#include <halmd/mdsim/gpu/forces/pair_short_ranged_kernel.cuh>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/utility/gpu/variant.cuh>

namespace halmd
{
namespace mdsim { namespace gpu { namespace forces
{
namespace lj_kernel
{

/** array of Lennard-Jones potential parameters for all combinations of particle types */
static texture<float4> param_;

/**
 * Lennard-Jones interaction of a pair of particles.
 */
class lj_potential
{
public:
    /**
     * Construct Lennard-Jones pair interaction potential.
     *
     * Fetch potential parameters from texture cache for particle pair.
     *
     * @param type1 type of first interacting particle
     * @param type2 type of second interacting particle
     */
    HALMD_GPU_ENABLED lj_potential(unsigned int type1, unsigned int type2)
      : pair_(
            tex1Dfetch(param_, symmetric_matrix::lower_index(type1, type2))
        ) {}

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
     * @returns tuple of absolute unit force and potential
     */
    template <typename float_type>
    HALMD_GPU_ENABLED tuple<float_type, float_type> operator()(float_type rr) const
    {
        float_type rri = pair_[SIGMA2] / rr;
        float_type ri6 = rri * rri * rri;
        float_type fval = 48 * pair_[EPSILON] * rri * ri6 * (ri6 - 0.5f) / pair_[SIGMA2];
        float_type en_pot = 4 * pair_[EPSILON] * ri6 * (ri6 - 1) - pair_[EN_CUT];

        return make_tuple(fval, en_pot);
    }

private:
    /** potential parameters for particle pair */
    fixed_vector<float, 4> pair_;
};

} // namespace lj_kernel

template <int dimension>
lj_wrapper<dimension> const lj_wrapper<dimension>::kernel = {
    lj_kernel::param_
};

template class lj_wrapper<3>;
template class lj_wrapper<2>;

template class pair_short_ranged_wrapper<3, lj_kernel::lj_potential>;
template class pair_short_ranged_wrapper<2, lj_kernel::lj_potential>;

}}} // namespace mdsim::gpu::forces

} // namespace halmd
