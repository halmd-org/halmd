/* CUDA device symbol
 *
 * Copyright (C) 2007  Peter Colberg
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

#ifndef CUDA_SYMBOL_HPP
#define CUDA_SYMBOL_HPP

#include <boost/noncopyable.hpp>
#include <cuda_runtime.h>

#ifndef __CUDACC__
# include <cuda_wrapper/error.hpp>
#endif

namespace cuda
{

/**
 * CUDA device symbol constant
 */
template <typename T>
class symbol
  : boost::noncopyable
{
public:
    typedef T value_type;
    typedef size_t size_type;

public:
    /**
     * initialize device symbol constant
     */
    explicit symbol(value_type const& symbol)
      : ptr_(&symbol)
    {}

    /**
     * return element count of device symbol
     */
    size_type size() const
    {
        return 1;
    }

    /**
     * returns device pointer to device symbol
     */
    value_type const* data() const
    {
        return ptr_;
    }

private:
    /** device symbol pointer */
    value_type const* ptr_;
};


/**
 * CUDA device symbol vector
 */
template <typename T>
class symbol<T[]>
  : boost::noncopyable
{
public:
    typedef T value_type;
    typedef size_t size_type;

public:
    /**
     * initialize device symbol vector
     */
    template <typename Array>
    explicit symbol(Array const& array)
      : ptr_(array)
      , size_(sizeof(array) / sizeof(value_type))
    {}

    /**
     * return element count of device symbol
     */
    size_type size() const
    {
        return size_;
    }

    /**
     * returns device pointer to device symbol
     */
    value_type const* data() const
    {
        return ptr_;
    }

private:
    /** device symbol pointer */
    value_type const* ptr_;
    /** array size */
    size_type size_;
};

} // namespace cuda

#endif /* ! CUDA_SYMBOL_HPP */
