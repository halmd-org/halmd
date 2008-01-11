/* cuda_wrapper/vector.hpp
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
#include <cuda_wrapper/allocator.hpp>
#include <cuda_wrapper/symbol.hpp>
#include <cuda_wrapper/host/vector.hpp>
#include <algorithm>
#include <assert.h>

namespace cuda
{

namespace host
{

template <typename T>
class vector;

}

template <typename T>
class symbol;


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

    vector(const host::vector<T>& src): n(0), ptr(NULL)
    {
	vector<T> dst(src.dim());
	dst = src;
	swap(*this, dst);
    }

    vector(const symbol<T> &src): n(0), ptr(NULL)
    {
	vector<T> dst(1);
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
	    // copy from device memory area to device memory area
	    CUDA_CALL(cudaMemcpy(ptr, v.get_ptr(), n * sizeof(T), cudaMemcpyDeviceToDevice));
	}
	return *this;
    }

    vector<T>& operator=(const host::vector<T>& v)
    {
	assert(v.dim() == n);
	// copy from host memory area to device memory area
	CUDA_CALL(cudaMemcpy(ptr, v.get_ptr(), n * sizeof(T), cudaMemcpyHostToDevice));
	return *this;
    }

    vector<T>& operator=(const symbol<T>& symbol)
    {
	assert(symbol.dim() == n);
	// copy from device symbol to device memory area
	CUDA_CALL(cudaMemcpyFromSymbol(ptr, reinterpret_cast<const char *>(symbol.get_ptr()), n * sizeof(T), 0, cudaMemcpyDeviceToDevice));
	return *this;
    }

    vector<T>& operator=(const T& value)
    {
	host::vector<T> v(n);
	*this = v = value;
	return *this;
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

#endif /* CUDA_ARRAY_HPP */
