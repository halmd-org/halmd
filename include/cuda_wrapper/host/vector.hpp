/* CUDA page-locked host memory vector
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

#ifndef CUDA_HOST_VECTOR_HPP
#define CUDA_HOST_VECTOR_HPP

#include <vector>

#include <cuda_wrapper/host/allocator.hpp>

namespace cuda {
namespace host {

/**
 * CUDA page-locked host memory vector
 */
template <typename T>
class vector : public std::vector<T, allocator<T> >
{
public:
    typedef allocator<T> _Alloc;
    typedef std::vector<T, allocator<T> > _Base;
    typedef typename _Base::value_type value_type;
    typedef typename _Base::iterator iterator;
    typedef typename _Base::const_iterator const_iterator;
    typedef typename _Base::reverse_iterator reverse_iterator;
    typedef typename _Base::const_reverse_iterator const_reverse_iterator;
    typedef typename _Base::reference reference;
    typedef typename _Base::const_reference const_reference;
    typedef typename _Base::pointer pointer;
    typedef typename _Base::difference_type difference_type;
    typedef typename _Base::size_type size_type;

public:
    /** creates an empty vector */
    vector(_Alloc const& alloc = _Alloc()) : _Base(alloc) {}
    /** creates a vector with n elements */
    vector(size_type n, T const& t = T(), _Alloc const& alloc = _Alloc()) : _Base(n, t, alloc) {}
    /** creates a vector with a copy of a range */
    template <class InputIterator>
    vector(InputIterator begin, InputIterator end, _Alloc const& alloc = _Alloc()) : _Base(begin, end, alloc) {}

#if (CUDART_VERSION >= 2020)
    /**
     * returns device pointer to page-locked host memory
     */
    operator value_type*()
    {
        void* ptr = NULL;
        CUDA_CALL(cudaHostGetDevicePointer(&ptr, _Base::data(), 0));
        return reinterpret_cast<value_type*>(ptr);
    }

    operator value_type const*() const
    {
        void* ptr = NULL;
        CUDA_CALL(cudaHostGetDevicePointer(&ptr, _Base::data(), 0));
        return reinterpret_cast<value_type const*>(ptr);
    }
#endif
};

} // namespace host
} // namespace cuda

#endif /* ! CUDA_HOST_VECTOR_HPP */
