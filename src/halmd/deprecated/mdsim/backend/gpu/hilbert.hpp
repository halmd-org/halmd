/* Hilbert spacing-filling curve kernel
 *
 * Copyright © 2008-2009  Peter Colberg
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

#ifndef HALMD_MDSIM_GPU_HILBERT_HPP
#define HALMD_MDSIM_GPU_HILBERT_HPP

#include <cuda_wrapper.hpp>

namespace halmd { namespace gpu
{

struct hilbert_base
{
    static cuda::symbol<float> box;
    static cuda::symbol<uint> depth;
};

template <int dimension>
struct hilbert;

template <>
struct hilbert<3> : public hilbert_base
{
    static cuda::function<void (float4 const*, uint*)> curve;
};

template <>
struct hilbert<2> : public hilbert_base
{
    static cuda::function<void (float4 const*, uint*)> curve;
};

}} // namespace halmd::gpu

#endif /* ! HALMD_MDSIM_GPU_HILBERT_HPP */
