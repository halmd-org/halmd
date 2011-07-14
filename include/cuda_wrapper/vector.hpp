/* CUDA global device memory vector
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

#ifndef CUDA_VECTOR_HPP
#define CUDA_VECTOR_HPP

#include <algorithm>
#include <boost/shared_ptr.hpp>
#include <boost/utility.hpp>

namespace cuda {

/**
 * CUDA global device memory vector
 */
template <typename T>
class vector
{
public:
    typedef vector<T> vector_type;
    typedef T value_type;
    typedef T* pointer;
    typedef T const* const_pointer;
    typedef size_t size_type;

private:
    class container : boost::noncopyable
    {
    public:
        /**
         * allocate global device memory
         */
        container(size_type size) : m_size(size), m_ptr(NULL)
        {
            // avoid zero-size allocation with cudaMalloc
            if (m_size) {
                void* ptr;
                CUDA_CALL(cudaMalloc(&ptr, m_size * sizeof(value_type)));
                m_ptr = reinterpret_cast<pointer>(ptr);
            }
        }

        /**
         * free global device memory
         */
        ~container() throw()
        {
            // avoid freeing NULL pointer with cudaFree
            if (m_size) {
                cudaFree(reinterpret_cast<void*>(m_ptr));
            }
        }

        size_type size() const
        {
            return m_size;
        }

        operator pointer()
        {
            return m_ptr;
        }

        operator const_pointer() const
        {
            return m_ptr;
        }

    private:
        size_type m_size;
        pointer m_ptr;
    };

public:
    /**
     * initialize device vector of given size
     */
    vector(size_type size = 0) : m_mem(new container(size)), m_size(size) {}

    /**
     * returns element count of device vector
     */
    size_type size() const
    {
        return m_size;
    }

    /**
     * returns capacity
     */
    size_type capacity() const
    {
        return m_mem->size();
    }

    /**
     * resize element count of device vector
     */
    void resize(size_type size)
    {
        this->reserve(size);
        m_size = size;
    }

    /**
     * allocate sufficient memory for specified number of elements
     */
    void reserve(size_type size)
    {
        if (size > m_mem->size()) {
            m_mem.reset();
            m_mem.reset(new container(size));
        }
    }

    /**
     * swap device memory with vector
     */
    void swap(vector_type& v)
    {
        m_mem.swap(v.m_mem);
        std::swap(m_size, v.m_size);
    }

    /**
     * returns device pointer to allocated device memory
     */
    pointer data()
    {
        return *m_mem;
    }

    operator pointer()
    {
        return *m_mem;
    }

    const_pointer data() const
    {
        return *m_mem;
    }

    operator const_pointer() const
    {
        return *m_mem;
    }

private:
    boost::shared_ptr<container> m_mem;
    size_type m_size;
};

} // namespace cuda

#endif /* CUDA_VECTOR_HPP */
