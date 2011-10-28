/*
 * Copyright © 2010-2011  Peter Colberg and Felix Höfling
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

// CUDA symbols (__constant__ or __device__ variables at module scope)
// in different translation units must not have the same name, therefore
// we need to compile all instantiations of mpc_kernel in the same
// translation unit.
//
// Although symbols are static to their translation unit, and symbols
// with the same name have different addresses, cudaMemcpyFromSymbol and
// cudaMemcpyToSymbol will *silently* fail to copy from or to the given
// symbol, as CUDA addresses a symbol by its mangled C++ name, not its
// address in memory.

#include <halmd/mdsim/gpu/potentials/lennard_jones_kernel.cuh>
#include <halmd/mdsim/gpu/potentials/lennard_jones_simple_kernel.cuh>
#include <halmd/mdsim/gpu/potentials/morse_kernel.cuh>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace potentials {

cuda::texture<float4> lennard_jones_wrapper::param = lennard_jones_kernel::param_;

cuda::symbol<float> lennard_jones_simple_wrapper::rr_cut = lennard_jones_simple_kernel::rr_cut_;
cuda::symbol<float> lennard_jones_simple_wrapper::en_cut = lennard_jones_simple_kernel::en_cut_;

cuda::texture<float4> morse_wrapper::param = morse_kernel::param_;
cuda::texture<float> morse_wrapper::rr_cut = morse_kernel::rr_cut_;

} // namespace potentials
} // namespace gpu
} // namespace mdsim
} // namespace halmd
