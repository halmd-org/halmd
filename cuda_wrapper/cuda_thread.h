/* cuda_thread.h
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

#ifndef __CUDA_THREAD_H__
#define __CUDA_THREAD_H__

#include "cuda_base.h"
#include "cuda_error.h"


/*
 * CUDA thread management
 */
class cuda_thread
{
public:
  /*
   * blocks until the device has completed all preceding requested tasks
   */
  static void synchronize()
  {
    CUDA_CALL(cudaThreadSynchronize());
  }

  /*
   * cleans up all runtime-related resources associated with calling thread
   */
  static void exit()
  {
    CUDA_CALL(cudaThreadExit());
  }
};


#endif /* ! __CUDA_THREAD_H__ */
