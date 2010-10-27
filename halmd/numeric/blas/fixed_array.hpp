/*
 * Copyright © 2008-2010  Peter Colberg
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

#ifndef HALMD_NUMERIC_BLAS_FIXED_ARRAY_HPP
#define HALMD_NUMERIC_BLAS_FIXED_ARRAY_HPP

#ifndef __CUDACC__
# include <boost/array.hpp>
#else
# include <cuda.h> // defines CUDA_VERSION
#endif

#if CUDA_VERSION <= 2030
# include <boost/type_traits/is_pod.hpp>
# include <boost/utility/enable_if.hpp>
# include <halmd/numeric/mp/dsfloat.hpp>
#endif

#include <halmd/config.hpp>

namespace halmd
{
namespace detail { namespace numeric { namespace blas
{

#ifndef __CUDACC__

template <typename T, size_t N>
struct fixed_array
  : boost::array<T, N> {};

#elif CUDA_VERSION <= 2030

template <typename T, typename Enable = void>
struct fixed_array_pod_type;

//
// The purpose of a fixed_array is to serve as the underlying
// array type to a fixed-length algebraic vector. It defines
// operator[] to allow convenient access of its components.
//
template <typename T, size_t N>
struct fixed_array
{
public:
    typedef T value_type;
    typedef value_type& reference;
    typedef value_type const& const_reference;
    typedef size_t size_type;
    enum { static_size = N };

    HALMD_GPU_ENABLED reference operator[](size_type i)
    {
        return (reinterpret_cast<value_type*>(this))[i];
    }

    HALMD_GPU_ENABLED const_reference operator[](size_type i) const
    {
        return (reinterpret_cast<value_type const*>(this))[i];
    }

private:
    typename fixed_array_pod_type<T>::type array[N];
};

//
// CUDA shared memory arrays only allow POD-type data members, so we
// cannot use structs with non-default constructors as a data members.
// Instead, we define an equivalent POD type here and cast to the
// non-POD type upon accessing an element of the fixed_array.
//

template <typename T>
struct fixed_array_pod_type<T,
  typename boost::enable_if<boost::is_pod<T> >::type>
{
    typedef T type;
};

template <>
struct fixed_array_pod_type<dsfloat>
{
    typedef float2 type;
};

#else /* __CUDACC__ */

//
// The purpose of a fixed_array is to serve as the underlying
// array type to a fixed-length algebraic vector. It defines
// operator[] to allow convenient access of its components.
//
// fixed_array models the Boost Collection concept.
//
template <typename T, size_t N>
class fixed_array
{
public:
    typedef T value_type;
    typedef value_type* iterator;
    typedef value_type const* const_iterator;
    typedef value_type& reference;
    typedef value_type const& const_reference;
    typedef value_type* pointer;
    typedef value_type const* const_pointer;
    typedef ssize_t difference_type;
    typedef size_t size_type;

    enum { static_size = N };

    HALMD_GPU_ENABLED reference operator[](size_type i)
    {
        return data_[i];
    }

    HALMD_GPU_ENABLED const_reference operator[](size_type i) const
    {
        return data_[i];
    }

    HALMD_GPU_ENABLED iterator begin()
    {
        return &((*this)[0]);
    }

    HALMD_GPU_ENABLED const_iterator begin() const
    {
        return &((*this)[0]);
    }

    HALMD_GPU_ENABLED iterator end()
    {
        return &((*this)[N - 1]);
    }

    HALMD_GPU_ENABLED const_iterator end() const
    {
        return &((*this)[N - 1]);
    }

    HALMD_GPU_ENABLED reference front()
    {
        return (*this)[0];
    }

    HALMD_GPU_ENABLED const_reference front() const
    {
        return (*this)[0];
    }

    HALMD_GPU_ENABLED reference back()
    {
        return (*this)[N - 1];
    }

    HALMD_GPU_ENABLED const_reference back() const
    {
        return (*this)[N - 1];
    }

    HALMD_GPU_ENABLED size_type size() const
    {
        return N;
    }

    HALMD_GPU_ENABLED bool empty() const
    {
        return N == 0;
    }

    HALMD_GPU_ENABLED void swap(fixed_array<T, N>& other)
    {
        fixed_array<T, N> temp = other;
        other = *this;
        *this = temp;
    }

private:
    value_type data_[N];
};

//
// fixed_array comparison operators
//

template <typename T, size_t N>
HALMD_GPU_ENABLED bool operator==(fixed_array<T, N> const& x, fixed_array<T, N> const& y)
{
    for (size_t i = 0; i < N; ++i) {
        if (!(x[i] == y[i])) {
            return false;
        }
    }
    return true;
}

template <typename T, size_t N>
HALMD_GPU_ENABLED bool operator!=(fixed_array<T, N> const& x, fixed_array<T, N> const& y)
{
    return !(x == y);
}

template <typename T, size_t N>
HALMD_GPU_ENABLED bool operator<(fixed_array<T, N> const& x, fixed_array<T, N> const& y)
{
    for (size_t i = 0; i < N; ++i) {
        if (!(x[i] < y[i])) {
            return false;
        }
    }
    return true;
}

template <typename T, size_t N>
HALMD_GPU_ENABLED bool operator>(fixed_array<T, N> const& x, fixed_array<T, N> const& y)
{
    return y < x;
}

template <typename T, size_t N>
HALMD_GPU_ENABLED bool operator<=(fixed_array<T, N> const& x, fixed_array<T, N> const& y)
{
    return !(y < x);
}

template <typename T, size_t N>
HALMD_GPU_ENABLED bool operator>=(fixed_array<T, N> const& x, fixed_array<T, N> const& y)
{
    return !(x < y);
}

#endif /* __CUDACC__ */

}}} // namespace detail::numeric::blas

} // namespace halmd

#endif /* ! HALMD_NUMERIC_BLAS_FIXED_ARRAY_HPP */
