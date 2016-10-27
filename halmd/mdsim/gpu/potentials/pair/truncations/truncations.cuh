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

#ifndef HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_TRUNCATIONS_CUH
#define HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_TRUNCATIONS_CUH

#include <boost/preprocessor/seq/for_each.hpp>

#include <halmd/mdsim/gpu/potentials/pair/truncations/force_shifted_kernel.cuh>
#include <halmd/mdsim/gpu/potentials/pair/truncations/sharp_kernel.cuh>
#include <halmd/mdsim/gpu/potentials/pair/truncations/shifted_kernel.cuh>
#include <halmd/mdsim/gpu/potentials/pair/truncations/smooth_r4_kernel.cuh>

#define HALMD_PAIR_POTENTIAL_TRUNCATIONS (smooth_r4)(sharp)(shifted)(force_shifted)

#define _HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_MAKE_WRAPPER(x) x##_wrapper
#define _HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_MAKE_KERNEL(x) x##_kernel

#define _HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_WRAPPER(r, kernel_type, truncation) \
    template class truncations::_HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_MAKE_WRAPPER(truncation)<kernel_type>;

#define HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_WRAPPERS(kernel_type) \
    BOOST_PP_SEQ_FOR_EACH(_HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_WRAPPER, kernel_type \
                        , HALMD_PAIR_POTENTIAL_TRUNCATIONS)

#define _HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_FORCE_KERNELS(r, kernel_type, truncation) \
    using namespace halmd::mdsim::gpu::potentials::pair::truncations::_HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_MAKE_KERNEL(truncation); \
    template class pair_trunc_wrapper<3, truncation<kernel_type> >; \
    template class pair_trunc_wrapper<2, truncation<kernel_type> >;


#define HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_FORCE_KERNELS(kernel_type) \
    BOOST_PP_SEQ_FOR_EACH(_HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_INSTANTIATE_FORCE_KERNELS, kernel_type\
                        , HALMD_PAIR_POTENTIAL_TRUNCATIONS)

#endif /* ! HALMD_MDSIM_GPU_POTENTIALS_PAIR_TRUNCATIONS_TRUNCATIONS_CUH */
