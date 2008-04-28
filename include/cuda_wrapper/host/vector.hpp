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

#ifndef CUDA_HOST_VECTOR_HPP
#define CUDA_HOST_VECTOR_HPP

#include <cuda/cuda_runtime.h>
#include <algorithm>
#include <vector>
#include <assert.h>
#include <cuda_wrapper/host/allocator.hpp>
#include <cuda_wrapper/vector.hpp>
#include <cuda_wrapper/symbol.hpp>
#include <cuda_wrapper/stream.hpp>
#include <cuda_wrapper/function.hpp>

namespace cuda
{

template <typename T>
class vector;

template <typename T>
class symbol;

namespace host
{

/**
 * vector container class for linear host memory
 */
template <typename T, typename Alloc = allocator<T> >
class vector : public std::vector<T, Alloc>
{
public:
    typedef Alloc _Alloc;
    typedef std::vector<T, Alloc> _Base;
    typedef vector<T, Alloc> vector_type;
    typedef T value_type;
    typedef size_t size_type;

public:
    /**
     * initialize host vector with copies of value
     */
    vector(size_type size, const value_type& value = value_type()) : _Base(size, value)
    {
    }

    /**
     * initialize host vector with copies of value
     */
    vector(config dim, const value_type& value = value_type()) : _Base(dim.threads(), value)
    {
    }

    /**
     * initialize host vector with content of host vector
     */
    vector(const vector_type& v) : _Base(v.size())
    {
	memcpy(v);
    }

    /**
     * initialize host vector with content of device vector
     */
    vector(const cuda::vector<value_type>& v) : _Base(v.size())
    {
	memcpy(v);
    }

    /**
     * copy from host memory area to host memory area
     */
    template <typename _Allocator>
    void memcpy(const vector<value_type, _Allocator>& v)
    {
	assert(v.size() == _Base::size());
	CUDA_CALL(cudaMemcpy(&_Base::front(), &v.front(), v.size() * sizeof(value_type), cudaMemcpyHostToHost));
    }

    /**
     * copy from device memory area to host memory area
     */
    void memcpy(const cuda::vector<value_type>& v)
    {
	assert(v.size() == _Base::size());
	CUDA_CALL(cudaMemcpy(&_Base::front(), v.data(), v.size() * sizeof(value_type), cudaMemcpyDeviceToHost));
    }

    /**
     * copy from device symbol to host memory area
     */
    void memcpy(const cuda::symbol<value_type>& symbol)
    {
	assert(symbol.size() == _Base::size());
	CUDA_CALL(cudaMemcpyFromSymbol(&_Base::front(), reinterpret_cast<const char *>(symbol.data()), symbol.size() * sizeof(value_type), 0, cudaMemcpyDeviceToHost));
    }

    /**
     * assign content of host vector to host vector
     */
    template <typename _Allocator>
    vector_type& operator=(const vector<value_type, _Allocator>& v)
    {
	if (this != &v) {
	    resize(v.size());
	    memcpy(v);
	}
	return *this;
    }

    /**
     * assign content of device vector to host vector
     */
    vector_type& operator=(const cuda::vector<value_type>& v)
    {
	resize(v.size());
	memcpy(v);
	return *this;
    }

    /**
     * assign content of device symbol to host vector
     */
    vector_type& operator=(const cuda::symbol<value_type>& symbol)
    {
	resize(symbol.size());
	memcpy(symbol);
	return *this;
    }

    /**
     * assign copies of value to host vector
     */
    vector_type& operator=(const value_type& value)
    {
	for (size_type i = 0; i < _Base::size(); i++) {
	    (*this)[i] = value;
	}
	return *this;
    }
};


/**
 * vector container class for linear page-locked host memory
 */
template <typename T>
class vector<T, allocator<T> > : public std::vector<T, allocator<T> >
{
public:
    typedef allocator<T> _Alloc;
    typedef std::vector<T, allocator<T> > _Base;
    typedef vector<T, allocator<T> > vector_type;
    typedef T value_type;
    typedef size_t size_type;

public:
    /**
     * initialize host vector with copies of value
     */
    vector(size_type size, const value_type& value = value_type()) : _Base(size, value)
    {
    }

    /**
     * initialize host vector with content of host vector
     */
    vector(const vector_type& v) : _Base(v.size())
    {
	memcpy(v);
    }

    /**
     * initialize host vector with content of device vector
     */
    vector(const cuda::vector<value_type>& v) : _Base(v.size())
    {
	memcpy(v);
    }

    /**
     * copy from host memory area to host memory area
     */
    template <typename _Allocator>
    void memcpy(const vector<value_type, _Allocator>& v)
    {
	assert(v.size() == _Base::size());
	CUDA_CALL(cudaMemcpy(&_Base::front(), &v.front(), v.size() * sizeof(value_type), cudaMemcpyHostToHost));
    }

    /**
     * copy from device memory area to host memory area
     */
    void memcpy(const cuda::vector<value_type>& v)
    {
	assert(v.size() == _Base::size());
	CUDA_CALL(cudaMemcpy(&_Base::front(), v.data(), v.size() * sizeof(value_type), cudaMemcpyDeviceToHost));
    }

    /**
     * copy from device symbol to host memory area
     */
    void memcpy(const cuda::symbol<value_type>& symbol)
    {
	assert(symbol.size() == _Base::size());
	CUDA_CALL(cudaMemcpyFromSymbol(&_Base::front(), reinterpret_cast<const char *>(symbol.data()), symbol.size() * sizeof(value_type), 0, cudaMemcpyDeviceToHost));
    }

#ifdef CUDA_WRAPPER_ASYNC_API

    /**
     * asynchronous copy from host memory area to host memory area
     *
     * requires page-locked host memory
     */
    void memcpy(const vector_type& v, const stream& stream)
    {
	assert(v.size() == _Base::size());
	CUDA_CALL(cudaMemcpyAsync(&_Base::front(), &v.front(), v.size() * sizeof(value_type), cudaMemcpyHostToHost, stream._stream));
    }

    /**
     * asynchronous copy from device memory area to host memory area
     *
     * requires page-locked host memory
     */
    void memcpy(const cuda::vector<value_type>& v, const stream& stream)
    {
	assert(v.size() == _Base::size());
	CUDA_CALL(cudaMemcpyAsync(&_Base::front(), v.data(), v.size() * sizeof(value_type), cudaMemcpyDeviceToHost, stream._stream));
    }

#endif /* CUDA_WRAPPER_ASYNC_API */

    /**
     * assign content of host vector to host vector
     */
    template <typename _Allocator>
    vector_type& operator=(const vector<value_type, _Allocator>& v)
    {
	if (this != &v) {
	    resize(v.size());
	    memcpy(v);
	}
	return *this;
    }

    /**
     * assign content of device vector to host vector
     */
    vector_type& operator=(const cuda::vector<value_type>& v)
    {
	resize(v.size());
	memcpy(v);
	return *this;
    }

    /**
     * assign content of device symbol to host vector
     */
    vector_type& operator=(const cuda::symbol<value_type>& symbol)
    {
	resize(symbol.size());
	memcpy(symbol);
	return *this;
    }

    /**
     * assign copies of value to host vector
     */
    vector_type& operator=(const value_type& value)
    {
	for (size_type i = 0; i < _Base::size(); i++) {
	    (*this)[i] = value;
	}
	return *this;
    }
};

}

}

#endif /* ! CUDA_HOST_VECTOR_HPP */
