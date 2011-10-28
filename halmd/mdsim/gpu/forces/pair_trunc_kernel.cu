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

// FIXME is this comment still relevant here?
//
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

#include <halmd/mdsim/gpu/forces/pair_trunc_kernel.cuh>

#include <halmd/mdsim/gpu/potentials/lennard_jones_kernel.cuh>
#include <halmd/mdsim/gpu/potentials/lennard_jones_simple_kernel.cuh>
#include <halmd/mdsim/gpu/potentials/morse_kernel.cuh>

namespace halmd {
namespace mdsim {
namespace gpu {
namespace forces {

using namespace potentials;

template class pair_trunc_wrapper<3, lennard_jones_kernel::lennard_jones>;
template class pair_trunc_wrapper<2, lennard_jones_kernel::lennard_jones>;

template class pair_trunc_wrapper<3, lennard_jones_simple_kernel::lennard_jones_simple>;
template class pair_trunc_wrapper<2, lennard_jones_simple_kernel::lennard_jones_simple>;

template class pair_trunc_wrapper<3, morse_kernel::morse>;
template class pair_trunc_wrapper<2, morse_kernel::morse>;

} // namespace forces
} // namespace gpu
} // namespace mdsim
} // namespace halmd
