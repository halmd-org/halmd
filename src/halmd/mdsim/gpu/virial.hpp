/* Parallel reduction kernel
 *
 * Copyright © 2008-2010  Peter Colberg
 *                        Felix Höfling
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

#ifndef HALMD_MDSIM_GPU_VIRIAL_HPP
#define HALMD_MDSIM_GPU_VIRIAL_HPP

#include <cuda_wrapper.hpp>
#include <halmd/math/gpu/dsfloat.cuh>
#include <halmd/math/gpu/dsvector.cuh>

namespace halmd { namespace gpu { namespace virial
{

enum {
    BLOCKS = 32,
    THREADS = 32 << DEVICE_SCALE
};

extern cuda::function<
    void(float4 const*, cu::vector<dsfloat, 4>*, uint),
    void(float2 const*, cu::vector<dsfloat, 2>*, uint),
    void(float4 const*, uint const*, cu::vector<dsfloat, 4>*, uint, uint),
    void(float2 const*, uint const*, cu::vector<dsfloat, 2>*, uint, uint)>
    sum;

}}} // namespace halmd::gpu::virial

#endif /* ! HALMD_MDSIM_GPU_VIRIAL_HPP */
