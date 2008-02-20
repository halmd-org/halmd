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
template <typename T, typename Alloc>
class vector;
#endif

}

#ifndef __CUDACC__
template <typename T>
class vector;
#endif


/**
 * vector container for a device symbol
 */
template <typename T>
class symbol
{
public:
    typedef symbol<T> vector_type;
    typedef T value_type;
    typedef size_t size_type;

protected:
    mutable size_type _size;
    const value_type* _ptr;

public:
    /**
     * initialize container with device symbol variable
     */
    symbol(const value_type& symbol) : _size(1), _ptr(&symbol)
    {
    }

    /**
     * initialize container with device symbol vector
     */
    symbol(const value_type* symbol) : _size(0), _ptr(symbol)
    {
    }

#ifndef __CUDACC__

    /**
     * copy from host memory area to device symbol
     */
    void memcpy(const host::vector<value_type>& v)
    {
	assert(v.size() == size());
	CUDA_CALL(cudaMemcpyToSymbol(reinterpret_cast<const char *>(get()), &v.front(), v.size() * sizeof(value_type), 0, cudaMemcpyHostToDevice));
    }

    /**
     * copy from device memory area to device symbol
     */
    void memcpy(const vector<value_type>& v)
    {
	assert(v.size() == size());
	CUDA_CALL(cudaMemcpyToSymbol(reinterpret_cast<const char *>(get()), v.get(), v.size() * sizeof(value_type), 0, cudaMemcpyDeviceToDevice));
    }

    /**
     * assign content of host vector to device symbol
     */
    vector_type& operator=(const host::vector<value_type>& v)
    {
	memcpy(v);
	return *this;
    }

    /**
     * assign content of device vector to device symbol
     */
    vector_type& operator=(const vector<value_type>& v)
    {
	memcpy(v);
	return *this;
    }

    /**
     * assign copies of value to device symbol
     */
    vector_type& operator=(const value_type& value)
    {
	host::vector<value_type> v(size(), value);
	memcpy(v);
	return *this;
    }

    /**
     * return element count of device symbol vector
     */
    size_type size() const
    {
	if (!_size) {
	    /*
	     * It would be preferable to issue the following CUDA runtime
	     * call directly upon construction. However, the constructor
	     * has to be compilable by the NVIDIA CUDA compiler as well,
	     * which does not support C++ runtime functionality, e.g.
	     * exceptions.
	     */
	    CUDA_CALL(cudaGetSymbolSize(&_size, reinterpret_cast<const char *>(get())));
	    _size /= sizeof(value_type);
	}

	return _size;
    }

#endif /* ! __CUDACC__ */

    /**
     * returns device pointer to device symbol
     */
    const value_type* get() const
    {
	return _ptr;
    }

private:
    // disable default copy constructor
    symbol(const vector_type&);
    // disable default assignment operator
    vector_type& operator=(const vector_type&);
};

}

#endif /* ! CUDA_SYMBOL_HPP */
