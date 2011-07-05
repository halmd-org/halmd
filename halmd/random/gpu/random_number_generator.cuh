/*
 * Copyright © 2010  Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef HALMD_RANDOM_GPU_RANDOM_KERNEL_CUH
#define HALMD_RANDOM_GPU_RANDOM_KERNEL_CUH

#include <halmd/random/gpu/rand48_kernel.cuh>
#include <halmd/utility/gpu/variant.cuh>

namespace halmd {
namespace random {
namespace gpu {

using utility::gpu::variant;
using utility::gpu::set;
using utility::gpu::get;

typedef variant<
    set<
        rand48_rng
    >
> random_number_generator;

} // namespace random
} // namespace gpu
} // namespace halmd

#endif /* ! HALMD_RANDOM_GPU_RANDOM_KERNEL_CUH */
