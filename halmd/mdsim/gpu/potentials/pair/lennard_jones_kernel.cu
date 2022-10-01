/*
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
 * Copyright © 2020       Jaslo Ziska
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
#include <halmd/mdsim/gpu/potentials/pair/adapters/hard_core_kernel.cuh>
#include <halmd/mdsim/gpu/potentials/pair/lennard_jones_kernel.hpp>
#include <halmd/mdsim/gpu/potentials/pair/truncations/truncations.cuh>
#include <halmd/numeric/blas/blas.hpp>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {
namespace pair {
namespace lennard_jones_kernel {

__device__ void lennard_jones::fetch_param(
    unsigned int type1, unsigned int type2
  , unsigned int ntype1, unsigned int ntype2
)
{
    pair_ = tex1Dfetch<float2>(t_param_, type1 * ntype2 + type2);
}

} // namespace lennard_jones_kernel

HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_WRAPPERS(lennard_jones_kernel::lennard_jones);

template class adapters::hard_core_wrapper<lennard_jones_kernel::lennard_jones>;
HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_WRAPPERS(
    adapters::hard_core_kernel::hard_core<lennard_jones_kernel::lennard_jones>
);

} // namespace pair
} // namespace potentials

// explicit instantiation of force kernels
namespace forces {

using namespace halmd::mdsim::gpu::potentials::pair::lennard_jones_kernel;
using namespace halmd::mdsim::gpu::potentials::pair::adapters::hard_core_kernel;

template class pair_full_wrapper<3, lennard_jones>;
template class pair_full_wrapper<2, lennard_jones>;
HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_FORCE_KERNELS(lennard_jones);

template class pair_full_wrapper<3, hard_core<lennard_jones> >;
template class pair_full_wrapper<2, hard_core<lennard_jones> >;
HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_FORCE_KERNELS(hard_core<lennard_jones>);

} // namespace forces

} // namespace gpu
} // namespace mdsim
} // namespace halmd
