/*
 * Copyright © 2010  Peter Colberg
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

#ifndef HALMD_CONFIG_HPP
#define HALMD_CONFIG_HPP

#ifdef __CUDACC__

#define HALMD_GPU_ENABLED __device__

#define HALMD_GPU_USING(__gpu__, __host__) using __gpu__

#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ >= 130
# define HALMD_GPU_DOUBLE_PRECISION
#endif

#else /* __CUDACC__ */

#define HALMD_GPU_ENABLED

#define HALMD_GPU_USING(__gpu__, __host__) using __host__

#endif /* ! __CUDACC__ */

#endif /* ! HALMD_CONFIG_HPP */
