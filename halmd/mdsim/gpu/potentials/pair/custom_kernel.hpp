/*
 * Copyright © 2023 Felix Höfling
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

#ifndef HALMD_MDSIM_GPU_POTENTIALS_PAIR_CUSTOM_KERNEL_HPP
#define HALMD_MDSIM_GPU_POTENTIALS_PAIR_CUSTOM_KERNEL_HPP

#include <halmd/numeric/blas/blas.hpp>
#include <halmd/utility/tuple.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace pair {
namespace custom_kernel {

/**
 * indices of potential parameters in float4 array
 */
enum {
    // FIXME rename PARAM[2-3] as in custom.hpp
    SIGMA     /**< interaction range parameter, in MD units */
  , PARAM2     /**< second potential parameter, in MD units */
  , PARAM3     /**< third potential parameter, in MD units */
};

template <typename float_type>
HALMD_GPU_ENABLED static inline tuple<float_type, float_type> compute(
    float_type const& rr
  , float_type const& sigma
  , float_type const& param2
  , float_type const& param3
)
{
    // FIXME
    // put here the actual formulas for the potential energy (en_pot) and
    // the force divided by the pair distance (fval).
    // use float_type, sqrt(rr), sigma, etc.
    float_type fval = -sigma * param2;
    float_type en_pot = param3 * rr / 2;

    return make_tuple(fval, en_pot);
}

/**
 * custom potential for the interaction of a pair of particles.
 */
class custom
{
public:
    /**
     * Construct custom's pair interaction potential.
     */
    custom(cudaTextureObject_t t_param) : t_param_(t_param) {}

    /**
     * Fetch potential parameters from texture cache for particle pair.
     *
     * @param type1 type of first interacting particle
     * @param type2 type of second interacting particle
     */
    HALMD_GPU_ENABLED void fetch_param(
        unsigned int type1, unsigned int type2
      , unsigned int ntype1, unsigned int ntype2
    );

    /**
     * Compute force and potential for interaction.
     *
     * @param rr squared distance between particles
     * @returns tuple of unit "force" @f$ -U'(r)/r @f$ and potential @f$ U(r) @f$
     */
    template <typename float_type>
    HALMD_GPU_ENABLED tuple<float_type, float_type> operator()(float_type rr) const
    {
        return custom_kernel::compute(rr, param_[SIGMA], param_[PARAM2], param_[PARAM3]);
    }

private:
    /** potential parameters for particle pair */
    // NOTE: we can have max. 4 potential parameters, see float4 texture in get_gpu_potential()
    fixed_vector<float, 4> param_;
    cudaTextureObject_t t_param_;
};

} // namespace custom_kernel

struct custom_wrapper {};

} // namespace pair
} // namespace potentials
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_POTENTIALS_PAIR_CUSTOM_KERNEL_HPP */
