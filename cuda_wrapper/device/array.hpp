/* cuda_wrapper/device/array.hpp
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

#ifndef CUDA_ARRAY_HPP
#define CUDA_ARRAY_HPP

#include <cuda_runtime.h>
#include <cuda_wrapper/host/array.hpp>
#include <cuda_wrapper/device/symbol.hpp>
#include <algorithm>
#include <assert.h>

namespace cuda
{

namespace host
{

template <typename T>
class array;

}

namespace device
{

template <typename T>
class symbol;


template <typename T>
class array
{
protected:
    size_t n;
    T *ptr;

public:
    array(size_t n): n(n)
    {
	void *p;
	// allocate linear memory on the device
	CUDA_CALL(cudaMalloc(&p, n * sizeof(T)));
	ptr = reinterpret_cast<T *>(p);
    }

    array(const array<T>& src): n(0), ptr(NULL)
    {
	if (this != &src) {
	    array<T> dst(src.dim());
	    dst = src;
	    swap(*this, dst);
	}
    }

    array(const host::array<T>& src): n(0), ptr(NULL)
    {
	array<T> dst(src.dim());
	dst = src;
	swap(*this, dst);
    }

    array(const symbol<T> &src): n(0), ptr(NULL)
    {
	array<T> dst(1);
	dst = src;
	swap(*this, dst);
    }

    ~array()
    {
	if (ptr != NULL) {
	    // free device memory
	    CUDA_CALL(cudaFree(ptr));
	}
    }

    array<T>& operator=(const array<T>& array)
    {
	if (this != &array) {
	    assert(array.dim() == n);
	    // copy from device memory area to device memory area
	    CUDA_CALL(cudaMemcpy(ptr, array.get_ptr(), n * sizeof(T), cudaMemcpyDeviceToDevice));
	}
	return *this;
    }

    array<T>& operator=(const host::array<T>& array)
    {
	assert(array.dim() == n);
	// copy from host memory area to device memory area
	CUDA_CALL(cudaMemcpy(ptr, array.get_ptr(), n * sizeof(T), cudaMemcpyHostToDevice));
	return *this;
    }

    array<T>& operator=(const symbol<T>& symbol)
    {
	assert(symbol.dim() == n);
	// copy from device symbol to device memory area
	CUDA_CALL(cudaMemcpyFromSymbol(ptr, reinterpret_cast<const char *>(symbol.get_ptr()), n * sizeof(T), 0, cudaMemcpyDeviceToDevice));
	return *this;
    }

    array<T>& operator=(const T& value)
    {
	host::array<T> array(n);
	*this = array = value;
	return *this;
    }

    static void swap(array<T>& a, array<T>& b)
    {
	std::swap(a.n, b.n);
	std::swap(a.ptr, b.ptr);
    }

    size_t dim() const
    {
	return n;
    }

    T *get_ptr() const
    {
	return ptr;
    }
};

}

}

#endif /* CUDA_ARRAY_HPP */
