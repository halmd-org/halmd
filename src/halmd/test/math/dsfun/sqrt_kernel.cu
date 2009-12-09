/* DSFUN square-root test
 *
 * Copyright (C) 2009  Peter Colberg
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

#include "sqrt_kernel.hpp"

__global__ void __kernel_sqrt(dsfloat const* g_a, dsfloat* g_b)
{
    unsigned int gid = threadIdx.x + blockIdx.x * blockDim.x;
    g_b[gid] = sqrt(g_a[gid]);
}

cuda::function<void (dsfloat const*, dsfloat*)> kernel_sqrt(__kernel_sqrt);
