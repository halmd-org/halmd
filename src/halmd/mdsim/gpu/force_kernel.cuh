/*
 * Copyright © 2008-2010  Peter Colberg
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

#ifndef HALMD_MDSIM_GPU_FORCE_KERNEL_CUH
#define HALMD_MDSIM_GPU_FORCE_KERNEL_CUH

#include <halmd/numeric/blas/blas.hpp>

namespace halmd { namespace mdsim { namespace gpu { namespace force_kernel
{

/**
 * Trace and off-diagonal elements of distance tensor
 */
template <typename T>
__device__ inline fixed_vector<T, 4> virial_tensor(T rr, fixed_vector<T, 3> const& r)
{
    fixed_vector<T, 4> v;
    v[0] = rr;
    v[1] = r[1] * r[2];
    v[2] = r[2] * r[0];
    v[3] = r[0] * r[1];
    return v;
}

template <typename T>
__device__ inline fixed_vector<T, 2> virial_tensor(T rr, fixed_vector<T, 2> const& r)
{
    fixed_vector<T, 2> v;
    v[0] = rr;
    v[1] = r[0] * r[1];
    return v;
}

}}}} // namespace halmd::mdsim::gpu::force_kernel

#endif /* ! HALMD_MDSIM_GPU_FORCE_KERNEL_CUH */
