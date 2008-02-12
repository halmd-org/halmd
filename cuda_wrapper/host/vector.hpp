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

#include <cuda_runtime.h>
#include <algorithm>
#include <vector>
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

/**
 * vector container class for linear host memory
 */
template <typename T, typename Alloc = allocator<T> >
class vector : public std::vector<T, Alloc>
{
protected:
    typedef std::vector<T, Alloc > _Base;

public:
    /**
     * initialize host vector with copies of value
     */
    vector(size_t n, const T& value = T()) : _Base(n, value)
    {
    }

    /**
     * initialize host vector with content of host vector
     */
    vector(const vector<T>& v) : _Base(v.dim())
    {
	memcpy(v);
    }

    /**
     * initialize host vector with content of device vector
     */
    vector(const cuda::vector<T>& v) : _Base(v.dim())
    {
	memcpy(v);
    }

    /**
     * copy from host memory area to host memory area
     */
    void memcpy(const vector<T, Alloc>& v)
    {
	assert(v.size() == _Base::size());
	CUDA_CALL(cudaMemcpy(&_Base::front(), &v.front(), v.size() * sizeof(T), cudaMemcpyHostToHost));
    }

    /**
     * copy from device memory area to host memory area
     */
    void memcpy(const cuda::vector<T>& v)
    {
	assert(v.dim() == _Base::size());
	CUDA_CALL(cudaMemcpy(&_Base::front(), v.get_ptr(), v.dim() * sizeof(T), cudaMemcpyDeviceToHost));
    }

    /**
     * copy from device symbol to host memory area
     */
    void memcpy(const cuda::symbol<T>& symbol)
    {
	assert(symbol.dim() == _Base::size());
	CUDA_CALL(cudaMemcpyFromSymbol(&_Base::front(), reinterpret_cast<const char *>(symbol.get_ptr()), symbol.size() * sizeof(T), 0, cudaMemcpyDeviceToHost));
    }

    /**
     * assign content of host vector to host vector
     */
    vector<T, Alloc>& operator=(const vector<T, Alloc>& v)
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
    vector<T, Alloc>& operator=(const cuda::vector<T>& v)
    {
	resize(v.dim());
	memcpy(v);
	return *this;
    }

    /**
     * assign content of device symbol to host vector
     */
    vector<T, Alloc>& operator=(const cuda::symbol<T>& symbol)
    {
	resize(symbol.size());
	memcpy(symbol);
	return *this;
    }

    /**
     * assign copies of value to host vector
     */
    vector<T, Alloc>& operator=(const T& value)
    {
	for (size_t i = 0; i < _Base::size(); i++) {
	    (*this)[i] = value;
	}
	return *this;
    }
};

}

}

#endif /* ! CUDA_HOST_VECTOR_HPP */
