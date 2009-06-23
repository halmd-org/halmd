/* DSFUN division test
 *
 * Copyright (C) 2009  Peter Colberg
 *
 * This program is free software: you can redistribute it and/or modify
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

#include "div_kernel.hpp"

__global__ void __kernel_div(dsfloat const* g_a, dsfloat const* g_b, dsfloat* g_c)
{
    unsigned int gid = threadIdx.x + blockIdx.x * blockDim.x;
    g_c[gid] = g_a[gid] / g_b[gid];
}

cuda::function<void (dsfloat const*, dsfloat const*, dsfloat*)> kernel_div(__kernel_div);
