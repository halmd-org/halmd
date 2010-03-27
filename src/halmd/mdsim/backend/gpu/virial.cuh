/* Lennard-Jones fluid kernel
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

#ifndef HALMD_MDSIM_GPU_VIRIAL_CUH
#define HALMD_MDSIM_GPU_VIRIAL_CUH

#include <halmd/math/gpu/vector2d.cuh>
#include <halmd/math/gpu/vector3d.cuh>
#include <halmd/math/gpu/vector4d.cuh>

namespace halmd { namespace cu { namespace virial
{

/**
 * trace and off-diagonal elements of distance tensor
 */
template <typename T>
__device__ vector<T, 4> tensor(T rr, vector<T, 3> const& r)
{
    return vector<T, 4>(rr, r.y * r.z, r.z * r.x, r.x * r.y);
}

template <typename T>
__device__ vector<T, 2> tensor(T rr, vector<T, 2> const& r)
{
    return vector<T, 2>(rr, r.x * r.y);
}

}}} // namespace halmd::cu::virial

#endif /* HALMD_MDSIM_GPU_VIRIAL_CUH */
