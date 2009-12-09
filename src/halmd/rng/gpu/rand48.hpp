/* Parallelized rand48 random number generator for CUDA
 *
 * Copyright © 2007-2009  Peter Colberg
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

#ifndef HALMD_RNG_GPU_RAND48_HPP
#define HALMD_RNG_GPU_RAND48_HPP

#include <cuda_wrapper.hpp>
#include <halmd/rng/gpu/uint48.cuh>

namespace halmd { namespace gpu { namespace rand48
{

extern cuda::function<void (uint48*)> leapfrog;
extern cuda::function<void (uint48 const*, uint48 const*, uint48*, uint48*, uint)> set;
extern cuda::function<void (ushort3*)> save;
extern cuda::function<void (uint48 const*, uint48 const*, uint48*, uint48*, ushort3)> restore;
extern cuda::function<void (float*, uint)> uniform;
extern cuda::function<void (uint*, uint)> get;

extern cuda::symbol<uint48> a;
extern cuda::symbol<uint48> c;
extern cuda::symbol<ushort3*> state;

}}} // namespace halmd::gpu::rand48

#endif /* ! HALMD_RNG_GPU_RAND48_HPP */
