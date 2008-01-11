/* cuda_wrapper/host/vector.hpp
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

#include <cuda_runtime.h>
#include <algorithm>
#include <assert.h>
#include <cuda_wrapper/host/allocator.hpp>
#include <cuda_wrapper/vector.hpp>
#include <cuda_wrapper/symbol.hpp>

namespace cuda
{

template <typename T>
class vector;

template <typename T>
class symbol;

namespace host
{

template <typename T>
class vector
{
protected:
    size_t n;
    T *ptr;

public:
    vector(size_t n): n(n), ptr(allocator<T>().allocate(n)) { }

    vector(const vector<T>& src): n(0), ptr(NULL)
    {
	if (this != &src) {
	    vector<T> dst(src.dim());
	    dst = src;
	    swap(*this, dst);
	}
    }

    vector(const cuda::vector<T>& src): n(0), ptr(NULL)
    {
	vector<T> dst(src.dim());
	dst = src;
	swap(*this, dst);
    }

    vector(const cuda::symbol<T>& src): n(0), ptr(NULL)
    {
	vector<T> dst(src.dim());
	dst = src;
	swap(*this, dst);
    }

    ~vector()
    {
	if (ptr != NULL) {
	    allocator<T>().deallocate(ptr, n);
	}
    }

    vector<T>& operator=(const vector<T>& v)
    {
	if (this != &v) {
	    assert(v.dim() == n);
	    // copy from host memory area to host memory area
	    CUDA_CALL(cudaMemcpy(ptr, v.get_ptr(), n * sizeof(T), cudaMemcpyHostToHost));
	}
	return *this;
    }

    vector<T>& operator=(const cuda::vector<T>& v)
    {
	assert(v.dim() == n);
	// copy from device memory area to host memory area
	CUDA_CALL(cudaMemcpy(ptr, v.get_ptr(), n * sizeof(T), cudaMemcpyDeviceToHost));
	return *this;
    }

    vector<T>& operator=(const cuda::symbol<T>& symbol)
    {
	assert(symbol.dim() == n);
	// copy from device symbol to host memory area
	CUDA_CALL(cudaMemcpyFromSymbol(ptr, reinterpret_cast<const char *>(symbol.get_ptr()), n * sizeof(T), 0, cudaMemcpyDeviceToHost));
	return *this;
    }

    vector<T>& operator=(const T& value)
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

    static void swap(vector<T>& a, vector<T>& b)
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
