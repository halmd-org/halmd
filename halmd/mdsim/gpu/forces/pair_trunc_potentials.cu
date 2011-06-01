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

#include <halmd/mdsim/gpu/forces/lennard_jones_kernel.cuh>
#include <halmd/mdsim/gpu/forces/lennard_jones_simple_kernel.cuh>
#include <halmd/mdsim/gpu/forces/morse_kernel.cuh>
