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
#include <halmd/mdsim/gpu/potentials/modified_lennard_jones_kernel.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/numeric/pow.hpp>  // std::pow is not a device function
#include <halmd/utility/tuple.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace modified_lennard_jones_kernel {

/** array of Lennard-Jones potential parameters for all combinations of particle types */
static texture<float4> param_;
/** squares of potential cutoff radius and energy shift for all combinations of particle types */
static texture<float2> rr_en_cut_;

/**
 * Lennard-Jones interaction of a pair of particles.
 */
class modified_lennard_jones
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
    HALMD_GPU_ENABLED modified_lennard_jones(unsigned int type1, unsigned int type2)
      : pair_(
            tex1Dfetch(param_, symmetric_matrix::lower_index(type1, type2))
        )
      , pair_rr_en_cut_(
            tex1Dfetch(rr_en_cut_, symmetric_matrix::lower_index(type1, type2))
        ) {}

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
     *   - U'(r) / r &=& 4 r^{-2} \epsilon (\sigma/r)^{n} \left[ m (\sigma/r)^{m-n} - n \right] \\
     *   U(r) &=& 4 \epsilon (\sigma/r)^{n} \left[ (\sigma/r)^{m-n} - 1 \right] \\
     *   r \partial_r r \partial_r U(r) &=& 4 \epsilon (\sigma/r)^{n} \left[ m^2 (\sigma/r)^{m-n} - n^2 \right]
     * @f}
     */
    template <typename float_type>
    HALMD_GPU_ENABLED tuple<float_type, float_type, float_type> operator()(float_type rr) const
    {
        float_type rri = pair_[SIGMA2] / rr;
        unsigned short m_2 = static_cast<unsigned short>(pair_[INDEX_M_2]);
        unsigned short n_2 = static_cast<unsigned short>(pair_[INDEX_N_2]);
        float_type rni = halmd::pow(rri, n_2);
        float_type rmni = (m_2 - n_2 == n_2) ? rni : halmd::pow(rri, m_2 - n_2);
        float_type eps_rni = pair_[EPSILON] * rni;
        float_type fval = 8 * rri * eps_rni * (m_2 * rmni - n_2) / pair_[SIGMA2];
        float_type en_pot = 4 * eps_rni * (rmni - 1) - pair_rr_en_cut_[1];
        float_type hvir = 16 * eps_rni * (m_2 * m_2 * rmni - n_2 * n_2);

        return make_tuple(fval, en_pot, hvir);
    }

private:
    /** potential parameters for particle pair */
    fixed_vector<float, 4> pair_;
    /** squared cutoff radius and energy shift for particle pair */
    fixed_vector<float, 2> pair_rr_en_cut_;
};

} // namespace modified_lennard_jones_kernel

cuda::texture<float4> modified_lennard_jones_wrapper::param = modified_lennard_jones_kernel::param_;
cuda::texture<float2> modified_lennard_jones_wrapper::rr_en_cut = modified_lennard_jones_kernel::rr_en_cut_;

} // namespace potentials

// explicit instantiation of force kernels
namespace forces {

using potentials::modified_lennard_jones_kernel::modified_lennard_jones;

template class pair_full_wrapper<3, modified_lennard_jones>;
template class pair_full_wrapper<2, modified_lennard_jones>;

template class pair_trunc_wrapper<3, modified_lennard_jones>;
template class pair_trunc_wrapper<2, modified_lennard_jones>;

} // namespace forces

} // namespace gpu
} // namespace mdsim
} // namespace halmd
