/*
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

#include <cuda_wrapper/cuda_wrapper.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>

using halmd::dsfloat;
using halmd::detail::numeric::mp::dsmulss;

/**
 * DSFUN addition test
 */
__global__ void __kernel_add(dsfloat const* g_a, dsfloat const* g_b, dsfloat* g_c)
{
    unsigned int gid = threadIdx.x + blockIdx.x * blockDim.x;
    g_c[gid] = g_a[gid] + g_b[gid];
}

/**
 * DSFUN subtraction test
 */
__global__ void __kernel_sub(dsfloat const* g_a, dsfloat const* g_b, dsfloat* g_c)
{
    unsigned int gid = threadIdx.x + blockIdx.x * blockDim.x;
    g_c[gid] = g_a[gid] - g_b[gid];
}

/**
 * DSFUN multiplication test
 */
__global__ void __kernel_mul(dsfloat const* g_a, dsfloat const* g_b, dsfloat* g_c)
{
    unsigned int gid = threadIdx.x + blockIdx.x * blockDim.x;
    g_c[gid] = g_a[gid] * g_b[gid];
}

/**
 * DSFUN single-single multiplication test
 */
__global__ void __kernel_mulss(float const* g_a, float const* g_b, dsfloat* g_c)
{
    unsigned int gid = threadIdx.x + blockIdx.x * blockDim.x;
    dsfloat c;
    dsmulss(c.hi, c.lo, g_a[gid], g_b[gid]);
    g_c[gid] = c;
}

/**
 * DSFUN division test
 */
__global__ void __kernel_div(dsfloat const* g_a, dsfloat const* g_b, dsfloat* g_c)
{
    unsigned int gid = threadIdx.x + blockIdx.x * blockDim.x;
    g_c[gid] = g_a[gid] / g_b[gid];
}

/**
 * DSFUN square-root test
 */
__global__ void __kernel_sqrt(dsfloat const* g_a, dsfloat* g_b)
{
    unsigned int gid = threadIdx.x + blockIdx.x * blockDim.x;
    g_b[gid] = sqrt(g_a[gid]);
}

cuda::function<void (dsfloat const*, dsfloat const*, dsfloat*)> kernel_add(__kernel_add);
cuda::function<void (dsfloat const*, dsfloat const*, dsfloat*)> kernel_sub(__kernel_sub);
cuda::function<void (dsfloat const*, dsfloat const*, dsfloat*)> kernel_mul(__kernel_mul);
cuda::function<void (float const*, float const*, dsfloat*)> kernel_mulss(__kernel_mulss);
cuda::function<void (dsfloat const*, dsfloat const*, dsfloat*)> kernel_div(__kernel_div);
cuda::function<void (dsfloat const*, dsfloat*)> kernel_sqrt(__kernel_sqrt);
