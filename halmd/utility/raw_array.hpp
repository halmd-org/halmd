/*
 * Copyright © 2013 Felix Höfling
 * Copyright © 2015 Nicolas Höft
 * Copyright © 2012 Peter Colberg
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

#ifndef HALMD_UTILITY_RAW_ARRAY_HPP
#define HALMD_UTILITY_RAW_ARRAY_HPP

#include <cassert>
#include <cstddef>
#include <cstring>
#include <type_traits>

namespace halmd {

/**
 * Uninitialised heap-allocated fixed-size array.
 *
 * This data structure implements a drop-in replacement for STL std::vector,
 * with a fixed number of array elements and without value initialisation.
 *
 * raw_array provides a random-access container with minimal allocation time.
 */
template <typename T>
class raw_array
{
    static_assert(std::is_standard_layout<T>(), "raw_array is compatible with standard layout types only");
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
    explicit raw_array(size_type size) : capacity_(size), size_(size), storage_(allocate(size)) {}

    /**
     * Deallocate array.
     */
    ~raw_array()
    {
        deallocate(storage_);
    }

    /**
     * Construct empty array.
     */
    raw_array() : capacity_(0), size_(0), storage_(nullptr) {}

    /**
     * Move constructor.
     */
    raw_array(raw_array&& array) : capacity_(array.capacity_), size_(array.size_), storage_(array.storage_)
    {
        array.capacity_ = 0;
        array.size_ = 0;
        array.storage_ = nullptr;
    }

    /**
     * Move assignment.
     */
    raw_array& operator=(raw_array&& array)
    {
        if (&array != this) {
            deallocate(storage_);
            capacity_ = array.capacity_;
            size_ = array.size_;
            storage_ = array.storage_;
            array.capacity_ = 0;
            array.size_ = 0;
            array.storage_ = nullptr;
        }
        return *this;
    }

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
        assert(i < size_);
        return storage_[i];
    }

    /**
     * Returns const reference to array element at given offset.
     */
    const_reference operator[](size_type i) const
    {
        assert(i < size_);
        return storage_[i];
    }

    /**
     * allocate sufficient memory for specified number of elements
     */
    void reserve(size_type size)
    {
        if(size > capacity_) {
            if (size_ > 0) {
                pointer tmp_storage = allocate(size);
                std::memcpy(tmp_storage, storage_, sizeof(value_type) * size_);
                deallocate(storage_);
                storage_ = tmp_storage;
            } else {
                deallocate(storage_);
                storage_ = allocate(size);
            }
            capacity_ = size;
        }
    }

    /**
     * resize element count, new elements will be uninitialised
     */
    void resize(size_type size)
    {
        reserve(size);
        size_ = size;
    }

    /** deleted implicit copy constructor */
    raw_array(raw_array const&) = delete;
    /** deleted implicit assignment operator */
    raw_array& operator=(raw_array const&) = delete;

private:
    /** allocate uninitialised storage */
    static pointer allocate(size_type size)
    {
        return static_cast<pointer>(::operator new(size * sizeof(value_type)));
    }

    /** deallocate storage */
    static void deallocate(pointer p)
    {
        ::operator delete(p);
    }

    /** number of array elements memory reserved for */
    size_type capacity_;
    /** number of array elements */
    size_type size_;
    /** uninitialised storage */
    pointer storage_;
};

template<typename T>
typename raw_array<T>::iterator begin(raw_array<T>& a)
{
    return a.begin();
}

template<typename T>
typename raw_array<T>::const_iterator begin(raw_array<T> const& a)
{
    return a.begin();
}

template<typename T>
typename raw_array<T>::const_iterator cbegin(raw_array<T> const& a)
{
    return a.begin();
}

template<typename T>
typename raw_array<T>::iterator end(raw_array<T>& a)
{
    return a.end();
}

template<typename T>
typename raw_array<T>::const_iterator end(raw_array<T> const& a)
{
    return a.end();
}

template<typename T>
typename raw_array<T>::const_iterator cend(raw_array<T> const& a)
{
    return a.end();
}

} // namespace halmd

#endif //! HALMD_UTILITY_RAW_ARRAY_HPP
