/* cuda_wrapper/host/array.hpp
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

#ifndef CUDA_HOST_ARRAY_HPP
#define CUDA_HOST_ARRAY_HPP

#include <cuda/cuda_runtime.h>
#include <algorithm>
#include <assert.h>
#include <cuda_wrapper/device/array.hpp>
#include <cuda_wrapper/device/symbol.hpp>

namespace cuda
{

namespace device
{

template <typename T>
class array;

template <typename T>
class symbol;

}

namespace host
{

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
	// allocate page-locked host memory accessible to the device
	CUDA_CALL(cudaMallocHost(&p, n * sizeof(T)));
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

    array(const device::array<T>& src): n(0), ptr(NULL)
    {
	array<T> dst(src.dim());
	dst = src;
	swap(*this, dst);
    }

    array(const device::symbol<T>& src): n(0), ptr(NULL)
    {
	array<T> dst(src.dim());
	dst = src;
	swap(*this, dst);
    }

    ~array()
    {
	if (ptr != NULL) {
	    // free page-locked host memory
	    CUDA_CALL(cudaFreeHost(ptr));
	}
    }

    array<T>& operator=(const array<T>& array)
    {
	if (this != &array) {
	    assert(array.dim() == n);
	    // copy from host memory area to host memory area
	    CUDA_CALL(cudaMemcpy(ptr, array.get_ptr(), n * sizeof(T), cudaMemcpyHostToHost));
	}
	return *this;
    }

    array<T>& operator=(const device::array<T>& array)
    {
	assert(array.dim() == n);
	// copy from device memory area to host memory area
	CUDA_CALL(cudaMemcpy(ptr, array.get_ptr(), n * sizeof(T), cudaMemcpyDeviceToHost));
	return *this;
    }

    array<T>& operator=(const device::symbol<T>& symbol)
    {
	assert(symbol.dim() == n);
	// copy from device symbol to host memory area
	CUDA_CALL(cudaMemcpyFromSymbol(ptr, reinterpret_cast<const char *>(symbol.get_ptr()), n * sizeof(T), 0, cudaMemcpyDeviceToHost));
	return *this;
    }

    array<T>& operator=(const T& value)
    {
	for (size_t i = 0; i < n; i++) {
	    ptr[i] = value;
	}
	return *this;
    }

    T& operator[](const size_t i) const
    {
	assert(i < n);
	return ptr[i];
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

#endif /* ! CUDA_HOST_ARRAY_HPP */
