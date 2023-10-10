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

#include <halmd/mdsim/gpu/forces/pair_full_kernel.cuh>
#include <halmd/mdsim/gpu/forces/pair_trunc_kernel.cuh>
#include <halmd/mdsim/gpu/potentials/pair/yukawa_kernel.hpp>
#include <halmd/mdsim/gpu/potentials/pair/truncations/truncations.cuh>
#include <halmd/numeric/blas/blas.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace pair {
namespace yukawa_kernel {

__device__ void yukawa::fetch_param(
    unsigned int type1, unsigned int type2
  , unsigned int ntype1, unsigned int ntype2
)
{
    param_ = tex1Dfetch<float2>(t_param_, type1 * ntype2 + type2);
}

} // namespace yukawa_kernel

HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_WRAPPERS(yukawa_kernel::yukawa);

} // namespace pair
} // namespace potentials

// explicit instantiation of force kernels
namespace forces {

using namespace halmd::mdsim::gpu::potentials::pair::yukawa_kernel;

template class pair_full_wrapper<3, yukawa>;
template class pair_full_wrapper<2, yukawa>;
HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_FORCE_KERNELS(yukawa);

} // namespace forces

} // namespace gpu
} // namespace mdsim
} // namespace halmd
