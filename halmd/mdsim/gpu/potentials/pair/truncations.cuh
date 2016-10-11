/*
 * Copyright © 2016 Daniel Kirchner
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

#ifndef HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_CUH
#define HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_CUH

#include <halmd/mdsim/gpu/potentials/pair/force_shifted_kernel.cuh>
#include <halmd/mdsim/gpu/potentials/pair/sharp_kernel.cuh>
#include <halmd/mdsim/gpu/potentials/pair/shifted_kernel.cuh>
#include <halmd/mdsim/gpu/potentials/pair/smooth_r4_kernel.cuh>

#define HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_WRAPPERS(kernel_type) \
    template class force_shifted_wrapper<kernel_type>; \
    template class sharp_wrapper<kernel_type>; \
    template class shifted_wrapper<kernel_type>; \
    template class smooth_r4_wrapper<kernel_type>;

#define HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_FORCE_KERNELS(kernel_type) \
    using namespace halmd::mdsim::gpu::potentials::pair::smooth_r4_kernel; \
    using namespace halmd::mdsim::gpu::potentials::pair::sharp_kernel; \
    using namespace halmd::mdsim::gpu::potentials::pair::shifted_kernel; \
    using namespace halmd::mdsim::gpu::potentials::pair::force_shifted_kernel; \
    template class pair_trunc_wrapper<3, smooth_r4<kernel_type> >; \
    template class pair_trunc_wrapper<2, smooth_r4<kernel_type> >; \
    template class pair_trunc_wrapper<3, sharp<kernel_type> >; \
    template class pair_trunc_wrapper<2, sharp<kernel_type> >; \
    template class pair_trunc_wrapper<3, shifted<kernel_type> >; \
    template class pair_trunc_wrapper<2, shifted<kernel_type> >; \
    template class pair_trunc_wrapper<3, force_shifted<kernel_type> >; \
    template class pair_trunc_wrapper<2, force_shifted<kernel_type> >;

#endif /* ! HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_CUH */
