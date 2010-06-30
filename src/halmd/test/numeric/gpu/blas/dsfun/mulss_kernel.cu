/* DSFUN single-single multiplication test
 *
 * Copyright © 2009  Peter Colberg
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

#include "mulss_kernel.hpp"

using namespace halmd::numeric::gpu::blas::detail;

__global__ void __kernel_mulss(float const* g_a, float const* g_b, dsfloat* g_c)
{
    unsigned int gid = threadIdx.x + blockIdx.x * blockDim.x;
    dsfloat c;
    dsmulss(c.hi, c.lo, g_a[gid], g_b[gid]);
    g_c[gid] = c;
}

cuda::function<void (float const*, float const*, dsfloat*)> kernel_mulss(__kernel_mulss);
