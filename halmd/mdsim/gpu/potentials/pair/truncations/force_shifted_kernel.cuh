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

#ifndef HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_FORCE_SHIFTED_KERNEL_CUH
#define HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_FORCE_SHIFTED_KERNEL_CUH

#include <halmd/mdsim/gpu/forces/pair_full_kernel.cuh>
#include <halmd/mdsim/gpu/forces/pair_trunc_kernel.cuh>
#include <halmd/mdsim/gpu/potentials/pair/truncations/force_shifted_kernel.hpp>
#include <halmd/numeric/blas/blas.hpp>
#include <halmd/utility/tuple.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace pair {
namespace truncations {
namespace force_shifted_kernel {

template <typename parent_kernel>
__device__ void force_shifted<parent_kernel>::fetch_param(
    unsigned int type1, unsigned int type2
  , unsigned int ntype1, unsigned int ntype2
)
{
    parent_kernel::fetch_param(type1, type2, ntype1, ntype2);
    pair_ = tex1Dfetch<float4>(t_param_, type1 * ntype2 + type2);
}

template <typename parent_kernel>
template <typename float_type>
__device__ tuple<float_type, float_type> force_shifted<parent_kernel>::operator()(float_type rr) const
{
    float_type f_abs, en_pot;
    float_type r = sqrtf(rr);
    tie(f_abs, en_pot) = parent_kernel::operator()(rr);
    f_abs -= pair_[FORCE_CUT] / r;
    en_pot = en_pot - pair_[EN_CUT] + (r - pair_[R_CUT]) * pair_[FORCE_CUT];
    return make_tuple(f_abs, en_pot);
}

} // namespace force_shifted_kernel
} // namespace truncations
} // namespace pair
} // namespace potentials
} // namespace gpu
} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_FORCE_SHIFTED_KERNEL_CUH */
