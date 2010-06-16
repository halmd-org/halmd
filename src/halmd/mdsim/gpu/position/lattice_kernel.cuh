/*
 * Copyright © 2008-2010  Peter Colberg and Felix Höfling
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

#ifndef HALMD_MDSIM_GPU_POSITION_LATTICE_KERNEL_CUH
#define HALMD_MDSIM_GPU_POSITION_LATTICE_KERNEL_CUH

#include <cuda_wrapper.hpp>

namespace halmd { namespace mdsim { namespace gpu { namespace position
{

template <int dimension>
struct lattice_wrapper
{
    static cuda::function<void (float4*, uint, float)> fcc;
    static cuda::function<void (float4*, uint, float)> sc;
};

}}}} // namespace halmd::mdsim::gpu::position

#endif /* ! HALMD_MDSIM_GPU_POSITION_LATTICE_KERNEL_CUH */
