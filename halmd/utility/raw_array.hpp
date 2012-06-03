/*
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
    explicit raw_array(size_type size) : size_(size), storage_(allocate(size)) {}

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
    raw_array() : size_(0), storage_(nullptr) {}

    /**
     * Move constructor.
     */
    raw_array(raw_array&& array) : size_(array.size_), storage_(array.storage_)
    {
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
            size_ = array.size_;
            storage_ = array.storage_;
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

    /** number of array elements */
    size_type size_;
    /** uninitialised storage */
    T* storage_;
};

} // namespace halmd

#endif //! HALMD_UTILITY_RAW_ARRAY_HPP
