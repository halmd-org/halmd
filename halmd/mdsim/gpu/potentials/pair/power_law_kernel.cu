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
#include <halmd/mdsim/gpu/potentials/pair/shifted_kernel.cuh>
#include <halmd/mdsim/gpu/potentials/pair/smooth_r4_kernel.cuh>
#include <halmd/mdsim/gpu/potentials/pair/power_law_kernel.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/numeric/pow.hpp>  // std::pow is not a device function
#include <halmd/utility/tuple.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace pair {
namespace power_law_kernel {

/** array of potential parameters for all combinations of particle types */
static texture<float4> param_;

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
    {}

    /**
     * Compute force and potential for interaction.
     *
     * @param rr squared distance between particles
     * @returns tuple of unit "force" @f$ -U'(r)/r @f$ and potential @f$ U(r) @f$
     *
     * @f{eqnarray*}{
     *   - U'(r) / r &=& n r^{-2} \epsilon (r/\sigma)^{-n} \\
     *   U(r) &=& \epsilon (r/\sigma)^{-n}
     * @f}
     */
    template <typename float_type>
    HALMD_GPU_ENABLED tuple<float_type, float_type> operator()(float_type rr) const
    {
        return compute(rr, pair_[SIGMA2], pair_[EPSILON], static_cast<unsigned short>(pair_[INDEX]));
    }

private:
    /** potential parameters for particle pair */
    fixed_vector<float, 4> pair_;
};

} // namespace power_law_kernel

cuda::texture<float4> power_law_wrapper::param = power_law_kernel::param_;
template class smooth_r4_wrapper<power_law_kernel::power_law>;
template class shifted_wrapper<power_law_kernel::power_law>;

} // namespace pair
} // namespace potentials

// explicit instantiation of force kernels
namespace forces {

using namespace halmd::mdsim::gpu::potentials::pair::power_law_kernel;
using namespace halmd::mdsim::gpu::potentials::pair::smooth_r4_kernel;
using namespace halmd::mdsim::gpu::potentials::pair::shifted_kernel;

template class pair_full_wrapper<3, power_law>;
template class pair_full_wrapper<2, power_law>;
template class pair_trunc_wrapper<3, smooth_r4<power_law> >;
template class pair_trunc_wrapper<2, smooth_r4<power_law> >;
template class pair_trunc_wrapper<3, shifted<power_law> >;
template class pair_trunc_wrapper<2, shifted<power_law> >;

} // namespace forces

} // namespace gpu
} // namespace mdsim
} // namespace halmd
