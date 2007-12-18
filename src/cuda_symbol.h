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
#include "cuda_array.h"


template <typename T>
class cuda_array;


template <typename T>
class cuda_symbol : public cuda_base
{
protected:
  const char *symbol;

public:
  cuda_symbol(T *symbol) : symbol(reinterpret_cast<const char *>(symbol))
  {
  }

  cuda_symbol& operator=(const T& value)
  {
    CUDA_CALL(cudaMemcpyToSymbol(symbol, &value, sizeof(T), 0, cudaMemcpyHostToDevice));
    return *this;
  }

  cuda_symbol& operator=(const cuda_array<T>& array)
  {
    CUDA_CALL(cudaMemcpyToSymbol(symbol, array.dev_ptr, array.n * sizeof(T), 0, cudaMemcpyDeviceToDevice));
    return *this;
  }
};


#endif /* ! __CUDA_SYMBOL_H__ */
