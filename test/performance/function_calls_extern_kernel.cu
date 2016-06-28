/*
 * Copyright © 2011  Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <test/performance/function_calls_extern_kernel.hpp>

/**
 * An empty function
 */
static __global__ void noop(double) {}

/**
 * Launch noop CUDA kernel
 */
void launch_noop_kernel(dim3 grid, dim3 block, double dummy)
{
    noop<<<grid, block>>>(dummy);
}

cuda::function<void (double)> noop_kernel(&noop);
