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

#include <vector>

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


#ifndef CUDA_KERNEL_MIN_ARG_SIZE
#define CUDA_KERNEL_MIN_ARG_SIZE 128
#endif


/*
 * CUDA execution control
 */
class cuda_kernel : public cuda_base
{
protected:
  /* function to execute on the device */
  const char *entry;
  /* argument stack */
  std::vector<char> arg_stack;

public:
  template <typename T>
  cuda_kernel(T *entry) : entry(reinterpret_cast<const char *>(entry))
  {
    arg_stack.reserve(CUDA_KERNEL_MIN_ARG_SIZE);
  }

  /*
   * push arbitrary argument onto argument stack
   */
  template <typename T>
  void setup_argument(T arg)
  {
    int offset = arg_stack.size();
    arg_stack.resize(offset + sizeof(T));
    memcpy(&arg_stack[offset], &arg, sizeof(T));
  }

  /*
   * push device memory pointer onto argument stack
   */
  template <typename T>
  void setup_argument(const cuda_array<T>& arg)
  {
    setup_argument(arg.dev_ptr);
  }

  /*
   * launch kernel with given execution parameters
   */
  void launch(const cuda_dim& dim, size_t shared_mem = 0)
  {
    CUDA_CALL(cudaConfigureCall(dim.grid, dim.block, shared_mem, 0));
    CUDA_CALL(cudaSetupArgument(&arg_stack[0], arg_stack.size(), 0));
    CUDA_CALL(cudaLaunch(entry));

    /* clear argument stack for next launch */
    arg_stack.clear();
  }
};


#endif /* ! __CUDA_KERNEL_H__ */
