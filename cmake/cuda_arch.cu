/*
 * Copyright © 2012  Peter Colberg
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

#define eval(arch) concat(arch)
#define concat(arch) __CUDA_ARCH__ ## arch ## __

/**
 * Generates nvcc error: identifier "__CUDA_ARCH__xxx__" is undefined
 */
static __global__ void cuda_arch()
{
    eval(__CUDA_ARCH__);
}

int main()
{
    cuda_arch<<<1, 1>>>();
    return 0;
}
