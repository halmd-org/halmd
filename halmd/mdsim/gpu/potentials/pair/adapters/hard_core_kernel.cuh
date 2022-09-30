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

#ifndef HALMD_MDSIM_GPU_POTENTIALS_PAIR_ADAPTERS_HARD_CORE_KERNEL_CUH
#define HALMD_MDSIM_GPU_POTENTIALS_PAIR_ADAPTERS_HARD_CORE_KERNEL_CUH

#include <halmd/mdsim/gpu/potentials/pair/adapters/hard_core_kernel.hpp>
#include <halmd/numeric/blas/blas.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace pair {
namespace adapters {
namespace hard_core_kernel {

template <typename parent_kernel>
__device__ void hard_core<parent_kernel>::fetch(
    unsigned int type1, unsigned int type2
  , unsigned int ntype1, unsigned int ntype2
)
{
    parent_kernel::fetch(type1, type2, ntype1, ntype2);
    r_core_ = tex1Dfetch<float>(t_param_, type1 * ntype2 + type2);
}

} // namespace hard_core_kernel
} // namespace adapters
} // namespace pair
} // namespace potentials
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_POTENTIALS_PAIR_ADAPTERS_HARD_CORE_KERNEL_CUH */
