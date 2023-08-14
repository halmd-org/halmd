/*
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

#ifndef HALMD_MDSIM_GPU_POTENTIALS_PAIR_ADAPTERS_HARD_CORE_KERNEL_HPP
#define HALMD_MDSIM_GPU_POTENTIALS_PAIR_ADAPTERS_HARD_CORE_KERNEL_HPP

#include <cuda_wrapper/cuda_wrapper.hpp>

#include <halmd/utility/tuple.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace pair {
namespace adapters {
namespace hard_core_kernel {

template <typename parent_kernel>
class hard_core
  : public parent_kernel
{
public:
    /**
     * Construct Hard Core Adapter.
     */
    hard_core(parent_kernel const& parent, cudaTextureObject_t t_param) :
        parent_kernel(parent), t_param_(t_param) {}

    /**
     * Fetch core parameter from texture cache for particle pair.
     *
     * @param type1 type of first interacting particle
     * @param type2 type of second interacting particle
     */
    HALMD_GPU_ENABLED void fetch(
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
        float_type r = sqrtf(rr);
        float_type r_s = r - r_core_;
        float_type f_abs, en_pot;
        tie(f_abs, en_pot) = parent_kernel::operator()(r_s * r_s);
        f_abs *= r_s / r;
        return make_tuple(f_abs, en_pot);
    }

private:
    /** core radius */
    float r_core_;
    cudaTextureObject_t t_param_;
};

} // namespace hard_core_kernel

template<typename parent_kernel>
struct hard_core_wrapper {};

} // namespace adapters
} // namespace pair
} // namespace potentials
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_POTENTIALS_PAIR_ADAPTERS_HARD_CORE_KERNEL_HPP */
