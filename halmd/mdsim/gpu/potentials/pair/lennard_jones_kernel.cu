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
#include <halmd/mdsim/gpu/potentials/pair/force_shifted_kernel.cuh>
#include <halmd/mdsim/gpu/potentials/pair/lennard_jones_kernel.hpp>
#include <halmd/mdsim/gpu/potentials/pair/sharp_kernel.cuh>
#include <halmd/mdsim/gpu/potentials/pair/shifted_kernel.cuh>
#include <halmd/mdsim/gpu/potentials/pair/smooth_r4_kernel.cuh>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/utility/tuple.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace pair {
namespace lennard_jones_kernel {

/** array of Lennard-Jones potential parameters for all combinations of particle types */
static texture<float2> param_;

/**
 * Lennard-Jones interaction of a pair of particles.
 */
class lennard_jones
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
    HALMD_GPU_ENABLED lennard_jones(
        unsigned int type1, unsigned int type2
      , unsigned int ntype1, unsigned int ntype2
    )
      : pair_(tex1Dfetch(param_, type1 * ntype2 + type2))
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
        return compute(rr, pair_[SIGMA2], pair_[EPSILON]);
    }

private:
    /** potential parameters for particle pair */
    fixed_vector<float, 2> pair_;
};

} // namespace lennard_jones_kernel

cuda::texture<float2> lennard_jones_wrapper::param = lennard_jones_kernel::param_;
template class smooth_r4_wrapper<lennard_jones_kernel::lennard_jones>;
template class sharp_wrapper<lennard_jones_kernel::lennard_jones>;
template class shifted_wrapper<lennard_jones_kernel::lennard_jones>;
template class force_shifted_wrapper<lennard_jones_kernel::lennard_jones>;

} // namespace pair
} // namespace potentials

// explicit instantiation of force kernels
namespace forces {

using namespace halmd::mdsim::gpu::potentials::pair::lennard_jones_kernel;
using namespace halmd::mdsim::gpu::potentials::pair::smooth_r4_kernel;
using namespace halmd::mdsim::gpu::potentials::pair::sharp_kernel;
using namespace halmd::mdsim::gpu::potentials::pair::shifted_kernel;
using namespace halmd::mdsim::gpu::potentials::pair::force_shifted_kernel;

template class pair_full_wrapper<3, lennard_jones>;
template class pair_full_wrapper<2, lennard_jones>;
template class pair_trunc_wrapper<3, smooth_r4<lennard_jones> >;
template class pair_trunc_wrapper<2, smooth_r4<lennard_jones> >;
template class pair_trunc_wrapper<3, sharp<lennard_jones> >;
template class pair_trunc_wrapper<2, sharp<lennard_jones> >;
template class pair_trunc_wrapper<3, shifted<lennard_jones> >;
template class pair_trunc_wrapper<2, shifted<lennard_jones> >;
template class pair_trunc_wrapper<3, force_shifted<lennard_jones> >;
template class pair_trunc_wrapper<2, force_shifted<lennard_jones> >;

} // namespace forces

} // namespace gpu
} // namespace mdsim
} // namespace halmd
