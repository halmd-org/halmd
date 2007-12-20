/* cuda.h
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

#ifndef __CUDA_DEVICE_H__
#define __CUDA_DEVICE_H__

#include "cuda_base.h"
#include "cuda_error.h"


/*
 * CUDA device management
 */
class cuda_device
{
public:
  /*
   * returns number of devices available for execution
   */
  static int count()
  {
    int count;
    CUDA_CALL(cudaGetDeviceCount(&count));
    return count;
  }

  /*
   * set device on which the active host thread executes device code
   */
  static void set(int dev)
  {
    CUDA_CALL(cudaSetDevice(dev));
  }

  /*
   * get device on which the active host thread executes device code
   */
  static int get()
  {
    int dev;
    CUDA_CALL(cudaGetDevice(&dev));
    return dev;
  }
};


#endif /* ! __CUDA_DEVICE_H__ */
