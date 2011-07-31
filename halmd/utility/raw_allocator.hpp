/*
 * Copyright © 2008-2011  Peter Colberg and Felix Höfling
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

/**
 * Allocator that wraps "C" malloc with empty construct() and destroy() methods
 *
 * The allocator is useful for huge STL containers of POD-like structs with
 * empty constructor. It avoids the element-wise initialisation upon resize()
 * or during construction of the container.
 *
 * This file is derived from ext/malloc_allocator.h of the
 * GNU Standard C++ Library, which wraps "C" malloc.
 */

#ifndef HALMD_UTILITY_RAW_ALLOCATOR_HPP
#define HALMD_UTILITY_RAW_ALLOCATOR_HPP

#include <bits/functexcept.h>
#include <cstdlib>
#include <new>

namespace halmd {
namespace detail {
namespace utility {

using std::size_t;
using std::ptrdiff_t;

template<typename _Tp>
class raw_allocator
{
public:
    typedef size_t     size_type;
    typedef ptrdiff_t  difference_type;
    typedef _Tp*       pointer;
    typedef const _Tp* const_pointer;
    typedef _Tp&       reference;
    typedef const _Tp& const_reference;
    typedef _Tp        value_type;

    template<typename _Tp1>
    struct rebind
    {
        typedef raw_allocator<_Tp1> other;
    };

    raw_allocator() throw() { }

    raw_allocator(const raw_allocator&) throw() { }

    template<typename _Tp1>
    raw_allocator(const raw_allocator<_Tp1>&) throw() { }

    ~raw_allocator() throw() { }

    pointer address(reference __x) const
    {
        return &__x;
    }

    const_pointer address(const_reference __x) const
    {
        return &__x;
    }

    // NB: __n is permitted to be 0.  The C++ standard says nothing
    // about what the return value is when __n == 0.
    pointer allocate(size_type __n, const void* = 0)
    {
        if (__builtin_expect(__n > this->max_size(), false))
            std::__throw_bad_alloc();

        pointer __ret = static_cast<_Tp*>(std::malloc(__n * sizeof(_Tp)));

        if (!__ret)
            std::__throw_bad_alloc();
        return __ret;
    }

    // __p is not permitted to be a null pointer.
    void deallocate(pointer __p, size_type) throw() // no-throw guarantee
    {
        std::free(static_cast<void*>(__p));
    }

    size_type max_size() const throw()
    {
        return size_t(-1) / sizeof(_Tp);
    }

    // memory at __p remains uninitialised
    void construct(pointer __p, const _Tp& __val) {}

    void destroy(pointer __p) {} // nothing to destroy
};

template<typename _Tp>
inline bool operator==(const raw_allocator<_Tp>&, const raw_allocator<_Tp>&)
{
    return true;
}

template<typename _Tp>
inline bool operator!=(const raw_allocator<_Tp>&, const raw_allocator<_Tp>&)
{
    return false;
}

} // namespace utility
} // namespace detail

// import into top-level namespace
using detail::utility::raw_allocator;

} // namespace halmd

#endif //! HALMD_UTILITY_RAW_ALLOCATOR_HPP
