/*
 * Copyright © 2013      Felix Höfling
 * Copyright © 2008-2010 Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#ifndef HALMD_NUMERIC_BLAS_DETAIL_SIZE_4_HPP
#define HALMD_NUMERIC_BLAS_DETAIL_SIZE_4_HPP

#include <halmd/config.hpp>

#include <algorithm>
#include <cassert>
#include <initializer_list>
#include <type_traits>

// CUDA vector types for host compiler
#ifdef HALMD_WITH_GPU
# include <cuda_runtime.h>
#endif

#include <halmd/numeric/blas/detail/array.hpp>
#include <halmd/numeric/blas/detail/vector.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>

namespace halmd {
namespace numeric {
namespace blas {
namespace detail {

/**
 * Four-dimensional single precision floating-point vector
 */
template <>
struct fixed_vector<float, 4>
  : fixed_array<float, 4>
{
    typedef fixed_array<float, 4> _Base;
    typedef _Base::value_type value_type;
    enum { static_size = _Base::static_size };

    HALMD_GPU_ENABLED fixed_vector() {}

    /**
     * Assign values from initializer list.
     */
    template <typename U>
    HALMD_GPU_ENABLED fixed_vector(std::initializer_list<U> const& v,
      typename std::enable_if<std::is_convertible<U, float>::value>::type* dummy = 0)
    {
        assert( v.size() == _Base::size() );
        std::copy(v.begin(), v.end(), _Base::begin());
    }

    /**
     * Initialization by scalar
     */
    HALMD_GPU_ENABLED fixed_vector(float const& s)
    {
        for (size_t i = 0; i < static_size; ++i) {
            (*this)[i] = s;
        }
    }

    /**
     * Explicit conversion from vector of convertible element type
     */
    template <typename U>
    HALMD_GPU_ENABLED explicit fixed_vector(fixed_vector<U, 4> const& v,
      typename std::enable_if<std::is_convertible<U, float>::value>::type* dummy = 0)
    {
        for (size_t i = 0; i < static_size; ++i) {
            (*this)[i] = static_cast<value_type>(v[i]);
        }
    }

#ifdef HALMD_WITH_GPU

    /**
     * Convert from CUDA vector type
     */
    HALMD_GPU_ENABLED fixed_vector(float4 const& v)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
        (*this)[2] = v.z;
        (*this)[3] = v.w;
    }

    /**
     * Convert to CUDA vector type
     */
    HALMD_GPU_ENABLED operator float4() const
    {
        float4 v;
        v.x = (*this)[0];
        v.y = (*this)[1];
        v.z = (*this)[2];
        v.w = (*this)[3];
        return v;
    }

#endif /* HALMD_WITH_GPU */
};

/**
 * Four-dimensional unsigned integer vector
 */
template <>
struct fixed_vector<unsigned int, 4>
  : fixed_array<unsigned int, 4>
{
    typedef fixed_array<unsigned int, 4> _Base;
    typedef _Base::value_type value_type;
    enum { static_size = _Base::static_size };

    HALMD_GPU_ENABLED fixed_vector() {}

    /**
     * Assign values from initializer list.
     */
    template <typename U>
    HALMD_GPU_ENABLED fixed_vector(std::initializer_list<U> const& v,
      typename std::enable_if<std::is_convertible<U, unsigned int>::value>::type* dummy = 0)
    {
        assert( v.size() == _Base::size() );
        std::copy(v.begin(), v.end(), _Base::begin());
    }

    /**
     * Initialization by scalar
     */
    HALMD_GPU_ENABLED fixed_vector(unsigned int const& s)
    {
        for (size_t i = 0; i < static_size; ++i) {
            (*this)[i] = s;
        }
    }

    /**
     * Explicit conversion from vector of convertible element type
     */
    template <typename U>
    HALMD_GPU_ENABLED explicit fixed_vector(fixed_vector<U, 4> const& v,
      typename std::enable_if<std::is_convertible<U, unsigned int>::value>::type* dummy = 0)
    {
        for (size_t i = 0; i < static_size; ++i) {
            (*this)[i] = static_cast<value_type>(v[i]);
        }
    }

#ifdef HALMD_WITH_GPU

    /**
     * Convert from CUDA vector type
     */
    HALMD_GPU_ENABLED fixed_vector(uint4 const& v)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
        (*this)[2] = v.z;
        (*this)[3] = v.w;
    }

    /**
     * Convert to CUDA vector type
     */
    HALMD_GPU_ENABLED operator uint4() const
    {
        uint4 v;
        v.x = (*this)[0];
        v.y = (*this)[1];
        v.z = (*this)[2];
        v.w = (*this)[3];
        return v;
    }

#endif /* HALMD_WITH_GPU */
};

/**
 * Four-dimensional integer vector
 */
template <>
struct fixed_vector<int, 4>
  : fixed_array<int, 4>
{
    typedef fixed_array<int, 4> _Base;
    typedef _Base::value_type value_type;
    enum { static_size = _Base::static_size };

    HALMD_GPU_ENABLED fixed_vector() {}

    /**
     * Assign values from initializer list.
     */
    template <typename U>
    HALMD_GPU_ENABLED fixed_vector(std::initializer_list<U> const& v,
      typename std::enable_if<std::is_convertible<U, int>::value>::type* dummy = 0)
    {
        assert( v.size() == _Base::size() );
        std::copy(v.begin(), v.end(), _Base::begin());
    }

    /**
     * Initialization by scalar
     */
    HALMD_GPU_ENABLED fixed_vector(int const& s)
    {
        for (size_t i = 0; i < static_size; ++i) {
            (*this)[i] = s;
        }
    }

    /**
     * Explicit conversion from vector of convertible element type
     */
    template <typename U>
    HALMD_GPU_ENABLED explicit fixed_vector(fixed_vector<U, 4> const& v,
      typename std::enable_if<std::is_convertible<U, int>::value>::type* dummy = 0)
    {
        for (size_t i = 0; i < static_size; ++i) {
            (*this)[i] = static_cast<value_type>(v[i]);
        }
    }

#ifdef HALMD_WITH_GPU

    /**
     * Convert from CUDA vector type
     */
    HALMD_GPU_ENABLED fixed_vector(int4 const& v)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
        (*this)[2] = v.z;
        (*this)[3] = v.w;
    }

    /**
     * Convert to CUDA vector type
     */
    HALMD_GPU_ENABLED operator int4() const
    {
        int4 v;
        v.x = (*this)[0];
        v.y = (*this)[1];
        v.z = (*this)[2];
        v.w = (*this)[3];
        return v;
    }

#endif /* HALMD_WITH_GPU */
};

/**
 * Four-dimensional double-single precision floating-point vector
 */
template <>
struct fixed_vector<dsfloat, 4> : fixed_array<dsfloat, 4>
{
    typedef fixed_array<dsfloat, 4> _Base;
    typedef _Base::value_type value_type;
    enum { static_size = _Base::static_size };

    HALMD_GPU_ENABLED fixed_vector() {}

    /**
     * Assign values from initializer list.
     */
    template <typename U>
    HALMD_GPU_ENABLED fixed_vector(std::initializer_list<U> const& v,
      typename std::enable_if<std::is_convertible<U, dsfloat>::value>::type* dummy = 0)
    {
        assert( v.size() == _Base::size() );
        std::copy(v.begin(), v.end(), _Base::begin());
    }

    /**
     * Initialization by scalar
     */
    template <typename U>
    HALMD_GPU_ENABLED fixed_vector(U const& s,
      typename std::enable_if<std::is_convertible<U, dsfloat>::value>::type* dummy = 0)
    {
        for (size_t i = 0; i < static_size; ++i) {
            (*this)[i] = s;
        }
    }

    /**
     * Implicit conversion from vector of convertible element type
     */
    template <typename U>
    HALMD_GPU_ENABLED fixed_vector(fixed_vector<U, 4> const& v,
      typename std::enable_if<std::is_convertible<U, dsfloat>::value>::type* dummy = 0)
    {
        for (size_t i = 0; i < static_size; ++i) {
            (*this)[i] = static_cast<value_type>(v[i]);
        }
    }

    HALMD_GPU_ENABLED fixed_vector(fixed_vector<float, 4> const& v, fixed_vector<float, 4> const& w)
    {
        for (size_t i = 0; i < static_size; ++i) {
            (*this)[i] = value_type(v[i], w[i]);
        }
    }
};

/**
 * Four-dimensional double precision floating-point vector
 */
template <>
struct fixed_vector<double, 4>
  : fixed_array<double, 4>
{
    typedef fixed_array<double, 4> _Base;
    typedef _Base::value_type value_type;
    enum { static_size = _Base::static_size };

    HALMD_GPU_ENABLED fixed_vector() {}

    /**
     * Assign values from initializer list.
     */
    template <typename U>
    HALMD_GPU_ENABLED fixed_vector(std::initializer_list<U> const& v,
      typename std::enable_if<std::is_convertible<U, double>::value>::type* dummy = 0)
    {
        assert( v.size() == _Base::size() );
        std::copy(v.begin(), v.end(), _Base::begin());
    }

    /**
     * Initialization by scalar
     */
    HALMD_GPU_ENABLED fixed_vector(double const& s)
    {
        for (size_t i = 0; i < static_size; ++i) {
            (*this)[i] = s;
        }
    }

    /**
     * Explicit conversion from vector of convertible element type
     */
    template <typename U>
    HALMD_GPU_ENABLED explicit fixed_vector(fixed_vector<U, 4> const& v,
      typename std::enable_if<std::is_convertible<U, double>::value>::type* dummy = 0)
    {
        for (size_t i = 0; i < static_size; ++i) {
            (*this)[i] = static_cast<value_type>(v[i]);
        }
    }

#ifdef HALMD_GPU_DOUBLE_PRECISION

    /**
     * Convert from CUDA vector type
     */
    HALMD_GPU_ENABLED fixed_vector(double4 const& v)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
        (*this)[2] = v.z;
        (*this)[3] = v.w;
    }

    /**
     * Convert to CUDA vector type
     */
    HALMD_GPU_ENABLED operator double4() const
    {
        double4 v;
        v.x = (*this)[0];
        v.y = (*this)[1];
        v.z = (*this)[2];
        v.w = (*this)[3];
        return v;
    }

#endif /* HALMD_GPU_DOUBLE_PRECISION */
};

} // namespace detail
} // namespace blas
} // namespace numeric
} // namespace halmd

#endif /* ! HALMD_NUMERIC_BLAS_DETAIL_SIZE_4_HPP */
