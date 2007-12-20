/* cuda_error.h
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

/*
 * CUDA runtime error checking
 */

#ifndef __CUDA_ERROR_H__
#define __CUDA_ERROR_H__

#include "cuda_base.h"


#define CUDA_ERROR(err) throw cuda_error(err)

#define CUDA_CALL(x) if (cudaSuccess != x) CUDA_ERROR(cudaGetLastError())


/*
 * CUDA error handling
 */
class cuda_error
{
public:
  /* CUDA error */
  const cudaError_t errno;

  cuda_error(cudaError_t _errno): errno(_errno)
  {
  }

  /*
   * returns a message string for the CUDA error
   */
  const char* what() const throw()
  {
    return cudaGetErrorString(errno);
  }
};


#endif /* ! __CUDA_ERROR_H__ */
