/* CUDA page-locked host memory vector
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

#include <cuda_wrapper/host/allocator.hpp>
#include <vector>

namespace cuda { namespace host
{

/**
 * CUDA page-locked host memory vector
 */
template <typename T>
class vector : public std::vector<T, allocator<T> >
{
public:
    typedef std::vector<T, allocator<T> > _Base;
    typedef typename _Base::size_type size_type;

public:
    /** creates an empty vector */
    vector() {}
    /** creates a vector with n elements */
    vector(size_type n) : _Base(n) {}
    /** creates a vector with n copies of t */
    vector(size_type n, T const& t) : _Base(n, t) {}
    /** creates a vector with a copy of a range */
    template <class InputIterator>
    vector(InputIterator begin, InputIterator end) : _Base(begin, end) {}
};

}} // namespace cuda::host

#endif /* ! CUDA_HOST_VECTOR_HPP */
