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

#ifndef CUDA_VECTOR_HPP
#define CUDA_VECTOR_HPP

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

template <typename T, typename Alloc>
class vector;

}

template <typename T>
class symbol;


/**
 * vector pseudo-container for linear global device memory
 */
template <typename T>
class vector
{
protected:
    size_t _size;
    T *_ptr;

public:
    /**
     * initialize device vector of given size
     */
    vector(size_t n): _size(n), _ptr(allocator<T>().allocate(n))
    {
    }

    /**
     * initialize device vector with content of device vector
     */
    vector(const vector<T>& src): _size(0), _ptr(NULL)
    {
	// ensure deallocation of device memory in case of an exception
	vector<T> dst(src.size());
	dst.memcpy(src);
	swap(*this, dst);
    }

    /**
     * initialize device vector with content of host vector
     */
    template <typename Alloc>
    vector(const host::vector<T, Alloc>& src): _size(0), _ptr(NULL)
    {
	// ensure deallocation of device memory in case of an exception
	vector<T> dst(src.size());
	dst.memcpy(src);
	swap(*this, dst);
    }

    /**
     * initialize device vector with content of device symbol
     */
    vector(const symbol<T> &src): _size(0), _ptr(NULL)
    {
	// ensure deallocation of device memory in case of an exception
	vector<T> dst(src.size());
	dst.memcpy(src);
	swap(*this, dst);
    }

    /**
     * deallocate device vector
     */
    ~vector()
    {
	if (_ptr != NULL) {
	    allocator<T>().deallocate(_ptr, _size);
	}
    }

    /**
     * copy from device memory area to device memory area
     */
    void memcpy(const vector<T>& v)
    {
	assert(v.size() == size());
	CUDA_CALL(cudaMemcpy(ptr(), v.ptr(), v.size() * sizeof(T), cudaMemcpyDeviceToDevice));
    }

    /*
     * copy from host memory area to device memory area
     */
    template <typename Alloc>
    void memcpy(const host::vector<T, Alloc>& v)
    {
	assert(v.size() == size());
	CUDA_CALL(cudaMemcpy(ptr(), &v.front(), v.size() * sizeof(T), cudaMemcpyHostToDevice));
    }

    /*
     * copy from device symbol to device memory area
     */
    void memcpy(const symbol<T>& symbol)
    {
	assert(symbol.size() == size());
	CUDA_CALL(cudaMemcpyFromSymbol(ptr(), reinterpret_cast<const char *>(symbol.ptr()), symbol.size() * sizeof(T), 0, cudaMemcpyDeviceToDevice));
    }

    /**
     * assign content of device vector to device vector
     */
    vector<T>& operator=(const vector<T>& v)
    {
	if (this != &v) {
	    memcpy(v);
	}
	return *this;
    }

    /**
     * assign content of host vector to device vector
     */
    template <typename Alloc>
    vector<T>& operator=(const host::vector<T, Alloc>& v)
    {
	memcpy(v);
	return *this;
    }

    /**
     * assign content of device symbol to device vector
     */
    vector<T>& operator=(const symbol<T>& symbol)
    {
	memcpy(symbol);
	return *this;
    }

    /**
     * assign copies of value to device vector
     */
    vector<T>& operator=(const T& value)
    {
	host::vector<T> v(size(), value);
	memcpy(v);
	return *this;
    }

    /**
     * swap device memory areas with another device vector
     */
    static void swap(vector<T>& a, vector<T>& b)
    {
	std::swap(a._size, b._size);
	std::swap(a._ptr, b._ptr);
    }

    /**
     * returns element count of device vector
     */
    size_t size() const
    {
	return _size;
    }

    /**
     * returns device pointer to allocated device memory
     */
    T *ptr() const
    {
	return _ptr;
    }
};

}

#endif /* CUDA_VECTOR_HPP */
