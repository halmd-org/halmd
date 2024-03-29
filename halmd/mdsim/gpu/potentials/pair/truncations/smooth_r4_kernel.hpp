/*
 * Copyright © 2012 Nicolas Höft
 * Copyright © 2016 Daniel Kirchner
 * Copyright © 2020 Jaslo Ziska
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

#ifndef HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_SMOOTH_R4_KERNEL_HPP
#define HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_SMOOTH_R4_KERNEL_HPP

#include <halmd/numeric/blas/blas.hpp>
#include <halmd/utility/tuple.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace pair {
namespace truncations {
namespace smooth_r4_kernel {

/**
 * indices of parameters
 */
enum {
    R_CUT       /**< cutoff distance */
  , RR_CUT      /**< square of cutoff distance */
  , EN_CUT      /**< potential energy at cutoff distance in MD units */
};

template <typename parent_kernel>
class smooth_r4
  : public parent_kernel
{
public:
    /**
     * Construct Smoothing Function.
     */
    smooth_r4(parent_kernel const& parent, cudaTextureObject_t t_param) :
        parent_kernel(parent), t_param_(t_param) {}

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
    HALMD_GPU_ENABLED tuple<float_type, float_type> operator()(float_type rr) const;

private:
    /** adapter parameters for particle pair */
    fixed_vector<float, 3> pair_;
    cudaTextureObject_t t_param_;
};

} // namespace smooth_r4_kernel

template<typename parent_kernel>
struct smooth_r4_wrapper
{
    static cuda::symbol<float> rri_smooth;
};

} // namespace truncations
} // namespace pair
} // namespace potentials
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_SMOOTH_R4_KERNEL_HPP */
