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

#include <cuda/cuda_runtime.h>
#include <cuda_wrapper/allocator.hpp>
#include <cuda_wrapper/symbol.hpp>
#include <cuda_wrapper/host/vector.hpp>
#include <cuda_wrapper/stream.hpp>
#include <cuda_wrapper/function.hpp>
#include <vector>
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
 * vector container for linear global device memory
 */
template <typename T>
class vector : protected std::vector<T, allocator<T> >
{
public:
    typedef allocator<T> _Alloc;
    typedef std::vector<T, allocator<T> > _Base;
    typedef vector<T> vector_type;
    typedef T value_type;
    typedef size_t size_type;

public:
    /**
     * initialize device vector of given size
     */
    vector(size_type size)
    {
	resize(size);
    }

    /**
     * initialize device vector of given size
     */
    vector(config dim)
    {
	resize(dim.threads());
    }

    /**
     * initialize device vector with content of device vector
     */
    vector(const vector_type& v)
    {
	resize(size);
	memcpy(v);
    }

    /**
     * initialize device vector with content of host vector
     */
    template <typename Alloc>
    vector(const host::vector<value_type, Alloc>& v)
    {
	resize(size);
	memcpy(v);
    }

    /**
     * initialize device vector with content of device symbol
     */
    vector(const symbol<value_type> &v)
    {
	resize(size);
	memcpy(v);
    }

    /**
     * copy from device memory area to device memory area
     */
    void memcpy(const vector_type& v)
    {
	assert(v.size() == size());
	CUDA_CALL(cudaMemcpy(data(), v.data(), v.size() * sizeof(value_type), cudaMemcpyDeviceToDevice));
    }

    /**
     * copy from host memory area to device memory area
     */
    template <typename Alloc>
    void memcpy(const host::vector<value_type, Alloc>& v)
    {
	assert(v.size() == size());
	CUDA_CALL(cudaMemcpy(data(), &v.front(), v.size() * sizeof(value_type), cudaMemcpyHostToDevice));
    }

#ifdef CUDA_WRAPPER_ASYNC_API

    /**
     * asynchronous copy from device memory area to device memory area
     */
    void memcpy(const vector_type& v, const stream& stream)
    {
	assert(v.size() == size());
	CUDA_CALL(cudaMemcpyAsync(data(), v.data(), v.size() * sizeof(value_type), cudaMemcpyDeviceToDevice, stream._stream));
    }

    /**
     * asynchronous copy from host memory area to device memory area
     *
     * requires page-locked host memory (default host vector allocator)
     */
    void memcpy(const host::vector<value_type, host::allocator<value_type> >& v, const stream& stream)
    {
	assert(v.size() == size());
	CUDA_CALL(cudaMemcpyAsync(data(), &v.front(), v.size() * sizeof(value_type), cudaMemcpyHostToDevice, stream._stream));
    }

#endif /* CUDA_WRAPPER_ASYNC_API */

    /*
     * copy from device symbol to device memory area
     */
    void memcpy(const symbol<value_type>& symbol)
    {
	assert(symbol.size() == size());
	CUDA_CALL(cudaMemcpyFromSymbol(data(), reinterpret_cast<const char *>(symbol.data()), symbol.size() * sizeof(value_type), 0, cudaMemcpyDeviceToDevice));
    }

    /**
     * assign content of device vector to device vector
     */
    vector_type& operator=(const vector_type& v)
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
    vector_type& operator=(const host::vector<value_type, Alloc>& v)
    {
	memcpy(v);
	return *this;
    }

    /**
     * assign content of device symbol to device vector
     */
    vector_type& operator=(const symbol<value_type>& symbol)
    {
	memcpy(symbol);
	return *this;
    }

    /**
     * assign copies of value to device vector
     */
    vector_type& operator=(const value_type& value)
    {
	host::vector<value_type, host::allocator<value_type> > v(size(), value);
	memcpy(v);
	return *this;
    }

    /**
     * returns element count of device vector
     */
    size_type size() const
    {
	return _Base::capacity();
    }

    /**
     * resize element count of device vector
     */
    void resize(size_type size)
    {
	// size of vector must always be kept at zero to prevent
	// initialization of vector elements, so use vector::reserve
	// in place of vector::resize
	_Base::reserve(size);
    }

    /**
     * returns device pointer to allocated device memory
     */
    value_type* data() const
    {
	return const_cast<value_type*>(&_Base::front());
    }
};

}

#endif /* CUDA_VECTOR_HPP */
