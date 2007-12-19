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

#include "union.h"

#include "cuda_base.h"
#include "cuda_error.h"
#include "cuda_array.h"


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
  template <typename X>
  static void setup_argument(const X& x)
  {
    CUDA_CALL(cudaSetupArgument(&x, sizeof(X), 0));
  }

  template <typename X, typename Y>
  static void setup_argument(const X& x, const Y& y)
  {
    union2<X, Y> args = { x, y };
    CUDA_CALL(cudaSetupArgument(&x, sizeof(X), __offsetof(args, x)));
    CUDA_CALL(cudaSetupArgument(&y, sizeof(Y), __offsetof(args, y)));
  }

  template <typename X, typename Y, typename Z>
  static void setup_argument(const X& x, const Y& y, const Z& z)
  {
    union3<X, Y, Z> args = { x, y, z };
    CUDA_CALL(cudaSetupArgument(&x, sizeof(X), __offsetof(args, x)));
    CUDA_CALL(cudaSetupArgument(&y, sizeof(Y), __offsetof(args, y)));
    CUDA_CALL(cudaSetupArgument(&z, sizeof(Z), __offsetof(args, z)));
  }

  template <typename X, typename Y, typename Z, typename W>
  static void setup_argument(const X& x, const Y& y, const Z& z, const W& w)
  {
    union4<X, Y, Z, W> args = { x, y, z, w };
    CUDA_CALL(cudaSetupArgument(&x, sizeof(X), __offsetof(args, x)));
    CUDA_CALL(cudaSetupArgument(&y, sizeof(Y), __offsetof(args, y)));
    CUDA_CALL(cudaSetupArgument(&z, sizeof(Z), __offsetof(args, z)));
    CUDA_CALL(cudaSetupArgument(&w, sizeof(W), __offsetof(args, w)));
  }

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
