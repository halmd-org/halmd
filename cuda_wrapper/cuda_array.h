/* cuda_array.h
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

#ifndef __CUDA_ARRAY_H__
#define __CUDA_ARRAY_H__

#include <assert.h>

#include "cuda_base.h"
#include "cuda_error.h"
#include "cuda_mem.h"
#include "cuda_host_array.h"

#ifndef __CUDACC__

template <typename T>
class cuda_host_array;

template <typename T>
class cuda_symbol;


template <typename T>
class cuda_array : public cuda_base
{
  friend class cuda_host_array<T>;
  friend class cuda_symbol<T>;

protected:
  T *dev_ptr;
  const int n;

public:
  cuda_array(int n): n(n)
  {
    cuda_mem::malloc(&dev_ptr, n);
  }

  ~cuda_array()
  {
    cuda_mem::free(dev_ptr);
  }

  cuda_array<T>& operator=(const cuda_array<T>& array)
  {
    assert(array.n == n);
    cuda_mem::DtoD(dev_ptr, array.dev_ptr, n);
    return *this;
  }

  cuda_array<T>& operator=(const cuda_host_array<T>& array)
  {
    assert(array.n == n);
    cuda_mem::HtoD(dev_ptr, array.host_ptr, n);
    return *this;
  }

  cuda_array<T>& operator=(const T& value)
  {
    cuda_host_array<T> array(n, value);
    return (*this = array);
  }

  int dim() const
  {
    return n;
  }

  T *front()
  {
    return dev_ptr;
  }

  T *at(int i)
  {
    return (dev_ptr + i);
  }
};

#endif /* ! __CUDACC__ */

#endif /* __CUDA_ARRAY_H__ */
