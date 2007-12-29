/* cuda_symbol.h
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

#ifndef __CUDA_SYMBOL_H__
#define __CUDA_SYMBOL_H__

#include "cuda_base.h"
#include "cuda_error.h"
#include "cuda_mem.h"
#include "cuda_array.h"
#include "cuda_host_array.h"


#ifndef __CUDACC__
template <typename T>
class cuda_array;
template <typename T>
class cuda_host_array;
#endif /* ! __CUDACC__ */


template <typename T>
class cuda_symbol : public cuda_base
{
protected:
  const T &symbol;
  const size_t n;

public:
  cuda_symbol(const T &symbol) : symbol(symbol), n(1) {}
  /* FIXME allow passing a symbol array */

#ifndef __CUDACC__
  cuda_symbol& operator=(const T& value)
  {
    assert(n == 1);
    cuda_mem::HtoS(symbol, &value, 1);
    return *this;
  }

  cuda_symbol& operator=(const cuda_host_array<T>& array)
  {
    assert(n == array.n);
    cuda_mem::HtoS(symbol, array.dev_ptr, n);
    return *this;
  }

  cuda_symbol& operator=(const cuda_array<T>& array)
  {
    assert(n == array.n);
    cuda_mem::DtoS(symbol, array.dev_ptr, n);
    return *this;
  }
#endif /* ! __CUDACC__ */
};


#endif /* ! __CUDA_SYMBOL_H__ */
