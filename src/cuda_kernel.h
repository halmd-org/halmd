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

#include <boost/preprocessor/repetition/enum_params.hpp>
#include <boost/preprocessor/repetition/enum_binary_params.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>
#include <boost/preprocessor/repetition/repeat_from_to.hpp>
#include <boost/preprocessor/inc.hpp>

#include "cuda_base.h"
#include "cuda_error.h"
#include "cuda_array.h"


/* maximum number of arguments passed to device functions */
#ifndef CUDA_KERNEL_MAX_ARGS
#define CUDA_KERNEL_MAX_ARGS 10
#endif


/*
 * CUDA execution configuration
 */
class cuda_dim
{
public:
  /* grid dimensions */
  const dim3 grid;
  /* block dimensions */
  const dim3 block;
  /* FIXME store useful numbers (no. of threads per grid/block) */

  cuda_dim(dim3 grid, dim3 block) : grid(grid), block(block)
  {
    /* FIXME store useful numbers (no. of threads per grid/block) */
  }

  int threads() const
  {
    return grid.y * grid.x * block.z * block.y * block.x;
  }

  int blocks_per_grid() const
  {
    return grid.y * grid.x;
  }

  int threads_per_block() const
  {
    return block.z * block.y * block.x;
  }
};


/*
 * CUDA execution control
 */
class cuda_kernel
{
public:
  /*
   * configure execution parameters
   */
  static void configure(const cuda_dim& dim, size_t shared_mem = 0)
  {
    CUDA_CALL(cudaConfigureCall(dim.grid, dim.block, shared_mem, 0));
  }

  /*
   * push arbitrary arguments into argument passing area
   */

#define DECL_STRUCT_MEMBER(z, n, text) \
  T##n x##n;

#define DECL_SETUP_ARGUMENT(z, n, text) \
  CUDA_CALL(cudaSetupArgument(&x##n, sizeof(T##n), offsetof(T_args, x##n)));

#define TEMPLATE(z, n, text) \
  template <BOOST_PP_ENUM_PARAMS_Z(z, n, typename T)> \
  static void setup_argument(BOOST_PP_ENUM_BINARY_PARAMS_Z(z, n, const T, &x)) \
  { \
    typedef struct { \
      BOOST_PP_REPEAT_##z(n, DECL_STRUCT_MEMBER, ) \
    } T_args; \
    BOOST_PP_REPEAT_##z(n, DECL_SETUP_ARGUMENT, ) \
  }

  BOOST_PP_REPEAT_FROM_TO(1, BOOST_PP_INC(CUDA_KERNEL_MAX_ARGS), TEMPLATE, )

#undef DECL_STRUCT_MEMBER
#undef DECL_SETUP_ARGUMENT
#undef TEMPLATE

  /*
   * launch kernel
   */
  template <typename T>
  static void launch(T *entry)
  {
    CUDA_CALL(cudaLaunch(reinterpret_cast<const char *>(entry)));
  }
};


#endif /* ! __CUDA_KERNEL_H__ */
