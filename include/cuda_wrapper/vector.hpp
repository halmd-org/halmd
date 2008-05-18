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
#include <cuda_wrapper/memory.hpp>
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
    vector() {}

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
    vector(const vector_type& src)
    {
	resize(size);
	copy(src, *this);
    }

    /**
     * initialize device vector with content of host vector
     */
    template <typename Alloc>
    vector(const host::vector<value_type, Alloc>& src)
    {
	resize(size);
	copy(src, *this);
    }

    /**
     * initialize device vector with content of device symbol
     */
    vector(const symbol<value_type> &src)
    {
	resize(size);
	copy(src, *this);
    }

    /**
     * assign content of device vector to device vector
     */
    vector_type& operator=(const vector_type& src)
    {
	if (this != &src) {
	    copy(src, *this);
	}
	return *this;
    }

    /**
     * assign content of host vector to device vector
     */
    template <typename Alloc>
    vector_type& operator=(const host::vector<value_type, Alloc>& src)
    {
	copy(src, *this);
	return *this;
    }

    /**
     * assign content of device symbol to device vector
     */
    vector_type& operator=(const symbol<value_type>& src)
    {
	copy(src, *this);
	return *this;
    }

    /**
     * assign copies of value to device vector
     */
    vector_type& operator=(const value_type& value)
    {
	host::vector<value_type, host::allocator<value_type> > src(size(), value);
	copy(src, *this);
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
    value_type* data()
    {
	return _Base::data();
    }

    /**
     * returns device pointer to allocated device memory
     */
    value_type const* data() const
    {
	return _Base::data();
    }
};


/**
 * returns device pointer to allocated device memory
 */
template <typename T>
T* cuda_cast(cuda::vector<T>& v)
{
    return v.data();
}

/**
 * returns device pointer to allocated device memory
 */
template <typename T>
T const* cuda_cast(cuda::vector<T> const& v)
{
    return v.data();
}

} // namespace cuda

#endif /* CUDA_VECTOR_HPP */
