/* cuda_kernel.h
 *
 * Copyright (C) 2007  Peter Colberg
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

#ifndef __CUDA_KERNEL_H__
#define __CUDA_KERNEL_H__

#include <boost/preprocessor/iteration/iterate.hpp>

/* maximum number of arguments passed to device functions */
#ifndef CUDA_KERNEL_MAX_ARGS
#define CUDA_KERNEL_MAX_ARGS 10
#endif

#define BOOST_PP_FILENAME_1 "cuda_kernel_template.h"
#define BOOST_PP_ITERATION_LIMITS (1, CUDA_KERNEL_MAX_ARGS)
#include BOOST_PP_ITERATE()

#endif /* ! __CUDA_KERNEL_H__ */
