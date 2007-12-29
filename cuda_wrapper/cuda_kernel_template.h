/* cuda_kernel_template.h
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

#include <boost/preprocessor/repetition/enum_params.hpp>
#include <boost/preprocessor/repetition/enum_binary_params.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>

#include "cuda_base.h"
#include "cuda_error.h"
#include "cuda_exec.h"

template <typename T>
struct cuda_kernel;

#define KERNEL_NUM_ARGS BOOST_PP_ITERATION()

#define KERNEL_TEMPLATE_PARMS BOOST_PP_ENUM_PARAMS(KERNEL_NUM_ARGS, typename T)

#define KERNEL_TEMPLATE_ARGS BOOST_PP_ENUM_PARAMS(KERNEL_NUM_ARGS, T)

#define KERNEL_PARMS BOOST_PP_ENUM_BINARY_PARAMS(KERNEL_NUM_ARGS, const T, &x)

#define KERNEL_SETUP_ARG(z, n, x) cuda_exec::setup_argument(x##n, &offset);


template <KERNEL_TEMPLATE_PARMS>
class cuda_kernel<void (KERNEL_TEMPLATE_ARGS)>
{
protected:
  typedef void T (KERNEL_TEMPLATE_ARGS);
  T *entry;

public:
  cuda_kernel(T *entry) : entry(entry) {}

#ifndef __CUDACC__
  void operator()(KERNEL_PARMS)
  {
    size_t offset = 0;
    BOOST_PP_REPEAT(KERNEL_NUM_ARGS, KERNEL_SETUP_ARG, x)
    cuda_exec::launch(entry);
  }
#endif /* ! __CUDACC__ */
};


#undef KERNEL_NUM_ARGS
#undef KERNEL_TEMPLATE_PARMS
#undef KERNEL_TEMPLATE_ARGS
#undef KERNEL_PARMS
#undef KERNEL_SETUP_ARG
