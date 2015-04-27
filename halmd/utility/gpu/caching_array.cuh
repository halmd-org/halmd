/*
 * Copyright © 2013,2015 Felix Höfling
 * Copyright © 2015      Nicolas Höft
 * Copyright © 2012      Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
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

#ifndef HALMD_UTILITLY_GPU_CACHING_ARRAY
#define HALMD_UTILITLY_GPU_CACHING_ARRAY

#include <boost/utility.hpp>
#include <cstddef>
#include <cub/util_allocator.cuh>
#include <cuda_wrapper/error.hpp>

namespace halmd {
namespace detail {

// the actual allocator for the device array, handling
// the allocations of the memory chunks and reusing allocated memory
extern cub::CachingDeviceAllocator caching_allocator_;

} // namespace detail

/**
 * Uninitialised cached GPU memory array.
 *
 * This data structure implements a drop-in replacement for STL std::vector
 * for GPU device memory using the CachingDeviceAllocator provided by CUB.
 *
 * The indended use is for allocation of (small) temporary memory chunks required
 * by device kernels. Note due to the restrictions by cub, this can be used in CUDA
 * source files (.cu) only.
 */
template <typename T>
class caching_array : boost::noncopyable
{
public:
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
    typedef T value_type;
    typedef T& reference;
    typedef T const& const_reference;
    typedef T* pointer;
    typedef T const* const_pointer;
    typedef T* iterator;
    typedef T const* const_iterator;

    /**
     * Allocate uninitialised array of given number of elements.
     */
    explicit caching_array(size_type size) : capacity_(size), size_(size), storage_(allocate(size)) {}

    /**
     * Deallocate array.
     */
    ~caching_array()
    {
        deallocate(storage_);
    }

    /**
     * Construct empty array.
     */
    caching_array() : capacity_(0), size_(0), storage_(0) {}

    /**
     * Returns iterator pointing to first element of array.
     */
    iterator begin()
    {
        return storage_;
    }

    /**
     * Returns const iterator pointing to first element of array.
     */
    const_iterator begin() const
    {
        return storage_;
    }

    /**
     * return capacity
     */
    size_type capacity() const
    {
        return capacity_;
    }

    /**
     * Set the array size to zero, does not affect the capcity of the array.
     */
    void clear()
    {
        size_ = 0;
    }

    /**
     * Returns iterator pointing one past last element of array.
     */
    iterator end()
    {
        return storage_ + size_;
    }

    /**
     * Returns const iterator pointing one past last element of array.
     */
    const_iterator end() const
    {
        return storage_ + size_;
    }

    /**
     * Returns number of array elements.
     */
    size_type size() const
    {
        return size_;
    }

    /**
     * Returns reference to array element at given offset.
     */
    reference operator[](size_type i)
    {
        return storage_[i];
    }

    /**
     * Returns const reference to array element at given offset.
     */
    const_reference operator[](size_type i) const
    {
        return storage_[i];
    }

    /**
     * Allocates sufficient memory for the specified number of elements,
     * copies the contents of the container up re-allocation of memory.
     */
    void reserve(size_type size)
    {
        deallocate(storage_);
        storage_ = allocate(size);
        capacity_ = size;
    }

    /**
     * resize element count, new elements will be uninitialised
     */
    void resize(size_type size)
    {
        reserve(size);
        size_ = size;
    }

    /**
     * Exchanges the contents of the container with those of 'other'.
     */
    void swap(caching_array& other)
    {
        std::swap(capacity_, other.capacity_);
        std::swap(size_,     other.size_);
        std::swap(storage_,  other.storage_);
    }

private:
    /** allocate uninitialised storage */
    static pointer allocate(size_type size)
    {
        void* p;
        CUDA_CALL(detail::caching_allocator_.DeviceAllocate(&p, size * sizeof(T)));
        return static_cast<pointer>(p);
    }

    /** deallocate storage */
    static void deallocate(pointer p)
    {
        CUDA_CALL(detail::caching_allocator_.DeviceFree(reinterpret_cast<void *>(p)));
    }

    /** number of array elements memory reserved for */
    size_type capacity_;
    /** number of array elements */
    size_type size_;
    /** uninitialised storage */
    pointer storage_;
};

template<typename T>
typename caching_array<T>::iterator begin(caching_array<T>& a)
{
    return a.begin();
}

template<typename T>
typename caching_array<T>::const_iterator begin(caching_array<T> const& a)
{
    return a.begin();
}

template<typename T>
typename caching_array<T>::const_iterator cbegin(caching_array<T> const& a)
{
    return a.begin();
}

template<typename T>
typename caching_array<T>::iterator end(caching_array<T>& a)
{
    return a.end();
}

template<typename T>
typename caching_array<T>::const_iterator end(caching_array<T> const& a)
{
    return a.end();
}

template<typename T>
typename caching_array<T>::const_iterator cend(caching_array<T> const& a)
{
    return a.end();
}

template<typename T>
void swap(caching_array<T>& lhs, caching_array<T>& rhs){
    lhs.swap(rhs);
}

} // namespace halmd

#endif /* ! HALMD_UTILITLY_GPU_CACHING_ARRAY */
