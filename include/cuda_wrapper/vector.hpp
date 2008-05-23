/* CUDA global device memory vector
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

#include <algorithm>
#include <cuda_wrapper/allocator.hpp>
#include <vector>


namespace cuda
{

/**
 * CUDA global device memory vector
 */
template <typename T>
class vector : private std::vector<T, allocator<T> >
{
public:
    typedef allocator<T> _Alloc;
    typedef std::vector<T, allocator<T> > _Base;
    typedef vector<T> vector_type;
    typedef T value_type;
    typedef size_t size_type;

public:
    vector() : size_(0) {}

    /**
     * initialize device vector of given size
     */
    vector(size_type size) : size_(size)
    {
	_Base::reserve(size_);
    }

    /**
     * returns element count of device vector
     */
    size_type size() const
    {
	return size_;
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
	size_ = size;
    }

    /**
     * returns capacity
     */
    size_type capacity() const
    {
	return _Base::capacity();
    }

    /**
     * allocate enough memory for specified number of elements
     */
    void reserve(size_type n)
    {
	_Base::reserve(std::max(size_, n));
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

private:
    // disable default copy constructor
    vector(vector_type const&);
    // disable default assignment operator
    vector_type& operator=(vector_type const&);

private:
    size_type size_;
};

} // namespace cuda

#endif /* CUDA_VECTOR_HPP */
