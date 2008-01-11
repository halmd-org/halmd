/* cuda_wrapper/symbol.hpp
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

#ifndef CUDA_SYMBOL_HPP
#define CUDA_SYMBOL_HPP

#include <cuda_runtime.h>
#ifndef __CUDACC__
#include <cuda_wrapper/error.hpp>
#include <cuda_wrapper/vector.hpp>
#include <cuda_wrapper/host/vector.hpp>
#endif

namespace cuda
{

namespace host
{

#ifndef __CUDACC__
template <typename T>
class vector;
#endif

}

#ifndef __CUDACC__
template <typename T>
class vector;
#endif


template <typename T>
class symbol
{
protected:
    const T *ptr;
    mutable size_t n;

public:
    /* constructor for device symbol variable */
    symbol(const T& symbol) : ptr(&symbol), n(1) {}

    /* constructor for device symbol vector */
    symbol(const T* symbol) : ptr(symbol), n(0) {}

#ifndef __CUDACC__
    symbol& operator=(const T& value)
    {
	host::vector<T> v(dim());
	*this = v = value;
	return *this;
    }

    symbol& operator=(const host::vector<T>& v)
    {
	assert(v.dim() == dim());
	// copy from host memory area to device symbol
	CUDA_CALL(cudaMemcpyToSymbol(reinterpret_cast<const char *>(ptr), v.get_ptr(), dim() * sizeof(T), 0, cudaMemcpyHostToDevice));
	return *this;
    }

    symbol& operator=(const vector<T>& v)
    {
	assert(v.dim() == dim());
	// copy from device memory area to device symbol
	CUDA_CALL(cudaMemcpyToSymbol(reinterpret_cast<const char *>(ptr), v.get_ptr(), dim() * sizeof(T), 0, cudaMemcpyDeviceToDevice));
	return *this;
    }

    size_t dim() const
    {
	if (n == 0) {
	    /*
	     * It would be preferable to issue the following CUDA runtime
	     * call directly upon construction. However, the constructor
	     * has to be compilable by the NVIDIA CUDA compiler as well,
	     * which does not support C++ runtime functionality, e.g.
	     * exceptions.
	     */
	    CUDA_CALL(cudaGetSymbolSize(&n, reinterpret_cast<const char *>(ptr)));
	    n = n / sizeof(T);
	}

	return n;
    }
#endif /* ! __CUDACC__ */

    T *get_ptr() const
    {
	return ptr;
    }
};

}

#endif /* ! CUDA_SYMBOL_HPP */
