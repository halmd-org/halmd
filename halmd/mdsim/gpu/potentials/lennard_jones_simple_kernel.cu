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

#include <halmd/algorithm/gpu/tuple.cuh>
#include <halmd/mdsim/gpu/forces/pair_full_kernel.cuh>
#include <halmd/mdsim/gpu/forces/pair_trunc_kernel.cuh>
#include <halmd/mdsim/gpu/potentials/lennard_jones_simple_kernel.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/utility/gpu/variant.cuh>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace lennard_jones_simple_kernel {

using algorithm::gpu::tuple;
using algorithm::gpu::make_tuple;

/** Lennard-Jones potential parameters: rr_cut, en_cut */
static __constant__ float rr_cut_;
static __constant__ float en_cut_;

/**
 * Lennard-Jones interaction for a simple fluid of a single species.
 */
class lennard_jones_simple
{
public:
    /**
     * Construct Lennard-Jones pair interaction potential.
     *
     * @param type1 type of first interacting particle
     * @param type2 type of second interacting particle
     */
    HALMD_GPU_ENABLED lennard_jones_simple(
        unsigned int type1, unsigned int type2
      , unsigned int ntype1, unsigned int ntype2
    )
    {}

    /**
     * Check whether particles are in interaction range.
     *
     * @param rr squared distance between particles
     */
    template <typename float_type>
    HALMD_GPU_ENABLED bool within_range(float_type rr) const
    {
        return (rr < rr_cut_);
    }

    /**
     * Compute force and potential for interaction.
     *
     * @param rr squared distance between particles
     * @returns tuple of unit "force" @f$ -U'(r)/r @f$, potential @f$ U(r) @f$,
     * and hypervirial @f$ r \partial_r r \partial_r U(r) @f$
     */
    template <typename float_type>
    HALMD_GPU_ENABLED tuple<float_type, float_type, float_type> operator()(float_type rr) const
    {
        const float_type sigma2 = 1;
        const float_type epsilon = 1;
        float_type rri = sigma2 / rr;
        float_type ri6 = rri * rri * rri;
        float_type eps_ri6 = epsilon * ri6;
        float_type fval = 48 * rri * eps_ri6 * (ri6 - 0.5f) / sigma2;
        float_type en_pot = 4 * eps_ri6 * (ri6 - 1) - en_cut_;
        float_type hvir = 576 * eps_ri6 * (ri6 - 0.25f);

        return make_tuple(fval, en_pot, hvir);
    }
};

} // namespace lennard_jones_simple_kernel

cuda::symbol<float> lennard_jones_simple_wrapper::rr_cut = lennard_jones_simple_kernel::rr_cut_;
cuda::symbol<float> lennard_jones_simple_wrapper::en_cut = lennard_jones_simple_kernel::en_cut_;

} // namespace potentials

// explicit instantiation of force kernels
namespace forces {

using potentials::lennard_jones_simple_kernel::lennard_jones_simple;

template class pair_full_wrapper<3, lennard_jones_simple>;
template class pair_full_wrapper<2, lennard_jones_simple>;

template class pair_trunc_wrapper<3, lennard_jones_simple>;
template class pair_trunc_wrapper<2, lennard_jones_simple>;

} // namespace forces

} // namespace gpu
} // namespace mdsim
} // namespace halmd
