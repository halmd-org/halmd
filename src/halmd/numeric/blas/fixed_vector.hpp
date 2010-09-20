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

#ifndef HALMD_NUMERIC_BLAS_FIXED_VECTOR_HPP
#define HALMD_NUMERIC_BLAS_FIXED_VECTOR_HPP

#include <boost/mpl/and.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/type_traits/is_floating_point.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/utility/enable_if.hpp>
#ifndef __CUDACC__
# include <cmath>
# include <iostream>
#endif
#ifdef WITH_CUDA
# include <cuda_runtime.h> // CUDA vector types for host compiler
#endif

#include <halmd/config.hpp>
#include <halmd/numeric/blas/fixed_array.hpp>

namespace halmd
{
namespace detail { namespace numeric { namespace blas
{

// import into current namespace
using namespace boost;

/**
 * N-dimensional vector of arbitrary value type
 */
template <typename T, size_t N>
struct fixed_vector
  : fixed_array<T, N>
{
    typedef fixed_array<T, N> _Base;
    typedef typename _Base::value_type value_type;
    enum { static_size = _Base::static_size };

    HALMD_GPU_ENABLED fixed_vector() {}

    /**
     * Initialization by scalar
     */
    template <typename U>
    HALMD_GPU_ENABLED fixed_vector(U const& s,
      typename enable_if<is_convertible<U, T> >::type* dummy = 0)
    {
        for (size_t i = 0; i < static_size; ++i) {
            (*this)[i] = s;
        }
    }

    /**
     * Explicit conversion from vector of convertible element type
     */
    template <typename U>
    HALMD_GPU_ENABLED explicit fixed_vector(fixed_vector<U, N> const& v,
      typename enable_if<is_convertible<U, T> >::type* dummy = 0)
    {
        for (size_t i = 0; i < static_size; ++i) {
            (*this)[i] = static_cast<value_type>(v[i]);
        }
    }

#ifndef __CUDACC__
    /**
     * implicit conversion from base classes
     */
    HALMD_GPU_ENABLED fixed_vector(fixed_array<T, N> const& v)
    {
        for (size_t i = 0; i < static_size; ++i) {
            (*this)[i] = v[i];
        }
    }
    HALMD_GPU_ENABLED fixed_vector(boost::array<T, N> const& v)
    {
        for (size_t i = 0; i < static_size; ++i) {
            (*this)[i] = v[i];
        }
    }
#endif
};

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
      typename enable_if<is_convertible<U, float> >::type* dummy = 0)
    {
        for (size_t i = 0; i < static_size; ++i) {
            (*this)[i] = static_cast<value_type>(v[i]);
        }
    }

#ifdef WITH_CUDA

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

#endif /* WITH_CUDA */
};

/**
 * Three-dimensional single precision floating-point vector
 */
template <>
struct fixed_vector<float, 3>
  : fixed_array<float, 3>
{
    typedef fixed_array<float, 3> _Base;
    typedef _Base::value_type value_type;
    enum { static_size = _Base::static_size };

    HALMD_GPU_ENABLED fixed_vector() {}

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
    HALMD_GPU_ENABLED explicit fixed_vector(fixed_vector<U, 3> const& v,
      typename enable_if<is_convertible<U, float> >::type* dummy = 0)
    {
        for (size_t i = 0; i < static_size; ++i) {
            (*this)[i] = static_cast<value_type>(v[i]);
        }
    }

#ifdef WITH_CUDA

    /**
     * Convert from CUDA vector type
     */
    HALMD_GPU_ENABLED fixed_vector(float3 const& v)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
        (*this)[2] = v.z;
    }

    HALMD_GPU_ENABLED fixed_vector(float4 const& v)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
        (*this)[2] = v.z;
    }

    /**
     * Convert to CUDA vector type
     */
    HALMD_GPU_ENABLED operator float3() const
    {
        float3 v;
        v.x = (*this)[0];
        v.y = (*this)[1];
        v.z = (*this)[2];
        return v;
    }

    HALMD_GPU_ENABLED operator float4() const
    {
        float4 v;
        v.x = (*this)[0];
        v.y = (*this)[1];
        v.z = (*this)[2];
        return v;
    }

#endif /* WITH_CUDA */
};

/**
 * Two-dimensional single precision floating-point vector
 */
template <>
struct fixed_vector<float, 2>
  : fixed_array<float, 2>
{
    typedef fixed_array<float, 2> _Base;
    typedef _Base::value_type value_type;
    enum { static_size = _Base::static_size };

    HALMD_GPU_ENABLED fixed_vector() {}

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
    HALMD_GPU_ENABLED explicit fixed_vector(fixed_vector<U, 2> const& v,
      typename enable_if<is_convertible<U, float> >::type* dummy = 0)
    {
        for (size_t i = 0; i < static_size; ++i) {
            (*this)[i] = static_cast<value_type>(v[i]);
        }
    }

#ifdef WITH_CUDA

    /**
     * Convert from CUDA vector type
     */
    HALMD_GPU_ENABLED fixed_vector(float2 const& v)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
    }

    HALMD_GPU_ENABLED fixed_vector(float3 const& v)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
    }

    HALMD_GPU_ENABLED fixed_vector(float4 const& v)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
    }

    /**
     * Convert to CUDA vector type
     */
    HALMD_GPU_ENABLED operator float2() const
    {
        float2 v;
        v.x = (*this)[0];
        v.y = (*this)[1];
        return v;
    }

    HALMD_GPU_ENABLED operator float3() const
    {
        float3 v;
        v.x = (*this)[0];
        v.y = (*this)[1];
        return v;
    }

    HALMD_GPU_ENABLED operator float4() const
    {
        float4 v;
        v.x = (*this)[0];
        v.y = (*this)[1];
        return v;
    }

#endif /* WITH_CUDA */
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
      typename enable_if<is_convertible<U, unsigned int> >::type* dummy = 0)
    {
        for (size_t i = 0; i < static_size; ++i) {
            (*this)[i] = static_cast<value_type>(v[i]);
        }
    }

#ifdef WITH_CUDA

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

#endif /* WITH_CUDA */
};

/**
 * Three-dimensional unsigned integer vector
 */
template <>
struct fixed_vector<unsigned int, 3>
  : fixed_array<unsigned int, 3>
{
    typedef fixed_array<unsigned int, 3> _Base;
    typedef _Base::value_type value_type;
    enum { static_size = _Base::static_size };

    HALMD_GPU_ENABLED fixed_vector() {}

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
    HALMD_GPU_ENABLED explicit fixed_vector(fixed_vector<U, 3> const& v,
      typename enable_if<is_convertible<U, unsigned int> >::type* dummy = 0)
    {
        for (size_t i = 0; i < static_size; ++i) {
            (*this)[i] = static_cast<value_type>(v[i]);
        }
    }

#ifdef WITH_CUDA

    /**
     * Convert from CUDA vector type
     */
    HALMD_GPU_ENABLED fixed_vector(uint3 const& v)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
        (*this)[2] = v.z;
    }

    HALMD_GPU_ENABLED fixed_vector(uint4 const& v)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
        (*this)[2] = v.z;
    }

    /**
     * Convert to CUDA vector type
     */
    HALMD_GPU_ENABLED operator uint3() const
    {
        uint3 v;
        v.x = (*this)[0];
        v.y = (*this)[1];
        v.z = (*this)[2];
        return v;
    }

    HALMD_GPU_ENABLED operator uint4() const
    {
        uint4 v;
        v.x = (*this)[0];
        v.y = (*this)[1];
        v.z = (*this)[2];
        return v;
    }

#endif /* WITH_CUDA */
};

/**
 * Two-dimensional unsigned integer vector
 */
template <>
struct fixed_vector<unsigned int, 2>
  : fixed_array<unsigned int, 2>
{
    typedef fixed_array<unsigned int, 2> _Base;
    typedef _Base::value_type value_type;
    enum { static_size = _Base::static_size };

    HALMD_GPU_ENABLED fixed_vector() {}

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
    HALMD_GPU_ENABLED explicit fixed_vector(fixed_vector<U, 2> const& v,
      typename enable_if<is_convertible<U, unsigned int> >::type* dummy = 0)
    {
        for (size_t i = 0; i < static_size; ++i) {
            (*this)[i] = static_cast<value_type>(v[i]);
        }
    }

#ifdef WITH_CUDA

    /**
     * Convert from CUDA vector type
     */
    HALMD_GPU_ENABLED fixed_vector(uint2 const& v)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
    }

    HALMD_GPU_ENABLED fixed_vector(uint3 const& v)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
    }

    HALMD_GPU_ENABLED fixed_vector(uint4 const& v)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
    }

    /**
     * Convert to CUDA vector type
     */
    HALMD_GPU_ENABLED operator uint2() const
    {
        uint2 v;
        v.x = (*this)[0];
        v.y = (*this)[1];
        return v;
    }

    HALMD_GPU_ENABLED operator uint3() const
    {
        uint3 v;
        v.x = (*this)[0];
        v.y = (*this)[1];
        return v;
    }

    HALMD_GPU_ENABLED operator uint4() const
    {
        uint4 v;
        v.x = (*this)[0];
        v.y = (*this)[1];
        return v;
    }

#endif /* WITH_CUDA */
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
      typename enable_if<is_convertible<U, int> >::type* dummy = 0)
    {
        for (size_t i = 0; i < static_size; ++i) {
            (*this)[i] = static_cast<value_type>(v[i]);
        }
    }

#ifdef WITH_CUDA

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

#endif /* WITH_CUDA */
};

/**
 * Three-dimensional integer vector
 */
template <>
struct fixed_vector<int, 3>
  : fixed_array<int, 3>
{
    typedef fixed_array<int, 3> _Base;
    typedef _Base::value_type value_type;
    enum { static_size = _Base::static_size };

    HALMD_GPU_ENABLED fixed_vector() {}

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
    HALMD_GPU_ENABLED explicit fixed_vector(fixed_vector<U, 3> const& v,
      typename enable_if<is_convertible<U, int> >::type* dummy = 0)
    {
        for (size_t i = 0; i < static_size; ++i) {
            (*this)[i] = static_cast<value_type>(v[i]);
        }
    }

#ifdef WITH_CUDA

    /**
     * Convert from CUDA vector type
     */
    HALMD_GPU_ENABLED fixed_vector(int3 const& v)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
        (*this)[2] = v.z;
    }

    HALMD_GPU_ENABLED fixed_vector(int4 const& v)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
        (*this)[2] = v.z;
    }

    /**
     * Convert to CUDA vector type
     */
    HALMD_GPU_ENABLED operator int3() const
    {
        int3 v;
        v.x = (*this)[0];
        v.y = (*this)[1];
        v.z = (*this)[2];
        return v;
    }

    HALMD_GPU_ENABLED operator int4() const
    {
        int4 v;
        v.x = (*this)[0];
        v.y = (*this)[1];
        v.z = (*this)[2];
        return v;
    }

#endif /* WITH_CUDA */
};

/**
 * Two-dimensional integer vector
 */
template <>
struct fixed_vector<int, 2>
  : fixed_array<int, 2>
{
    typedef fixed_array<int, 2> _Base;
    typedef _Base::value_type value_type;
    enum { static_size = _Base::static_size };

    HALMD_GPU_ENABLED fixed_vector() {}

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
    HALMD_GPU_ENABLED explicit fixed_vector(fixed_vector<U, 2> const& v,
      typename enable_if<is_convertible<U, int> >::type* dummy = 0)
    {
        for (size_t i = 0; i < static_size; ++i) {
            (*this)[i] = static_cast<value_type>(v[i]);
        }
    }

#ifdef WITH_CUDA

    /**
     * Convert from CUDA vector type
     */
    HALMD_GPU_ENABLED fixed_vector(int2 const& v)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
    }

    HALMD_GPU_ENABLED fixed_vector(int3 const& v)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
    }

    HALMD_GPU_ENABLED fixed_vector(int4 const& v)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
    }

    /**
     * Convert to CUDA vector type
     */
    HALMD_GPU_ENABLED operator int2() const
    {
        int2 v;
        v.x = (*this)[0];
        v.y = (*this)[1];
        return v;
    }

    HALMD_GPU_ENABLED operator int3() const
    {
        int3 v;
        v.x = (*this)[0];
        v.y = (*this)[1];
        return v;
    }

    HALMD_GPU_ENABLED operator int4() const
    {
        int4 v;
        v.x = (*this)[0];
        v.y = (*this)[1];
        return v;
    }

#endif /* WITH_CUDA */
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
     * Initialization by scalar
     */
    template <typename U>
    HALMD_GPU_ENABLED fixed_vector(U const& s,
      typename enable_if<is_convertible<U, dsfloat> >::type* dummy = 0)
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
      typename enable_if<is_convertible<U, dsfloat> >::type* dummy = 0)
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
 * Three-dimensional double-single precision floating-point vector
 */
template <>
struct fixed_vector<dsfloat, 3> : fixed_array<dsfloat, 3>
{
    typedef fixed_array<dsfloat, 3> _Base;
    typedef _Base::value_type value_type;
    enum { static_size = _Base::static_size };

    HALMD_GPU_ENABLED fixed_vector() {}

    /**
     * Initialization by scalar
     */
    template <typename U>
    HALMD_GPU_ENABLED fixed_vector(U const& s,
      typename enable_if<is_convertible<U, dsfloat> >::type* dummy = 0)
    {
        for (size_t i = 0; i < static_size; ++i) {
            (*this)[i] = s;
        }
    }

    /**
     * Implicit conversion from vector of convertible element type
     */
    template <typename U>
    HALMD_GPU_ENABLED fixed_vector(fixed_vector<U, 3> const& v,
      typename enable_if<is_convertible<U, dsfloat> >::type* dummy = 0)
    {
        for (size_t i = 0; i < static_size; ++i) {
            (*this)[i] = static_cast<value_type>(v[i]);
        }
    }

    HALMD_GPU_ENABLED fixed_vector(fixed_vector<float, 3> const& v, fixed_vector<float, 3> const& w)
    {
        for (size_t i = 0; i < static_size; ++i) {
            (*this)[i] = value_type(v[i], w[i]);
        }
    }
};

/**
 * Two-dimensional double-single precision floating-point vector
 */
template <>
struct fixed_vector<dsfloat, 2> : fixed_array<dsfloat, 2>
{
    typedef fixed_array<dsfloat, 2> _Base;
    typedef _Base::value_type value_type;
    enum { static_size = _Base::static_size };

    HALMD_GPU_ENABLED fixed_vector() {}

    /**
     * Initialization by scalar
     */
    template <typename U>
    HALMD_GPU_ENABLED fixed_vector(U const& s,
      typename enable_if<is_convertible<U, dsfloat> >::type* dummy = 0)
    {
        for (size_t i = 0; i < static_size; ++i) {
            (*this)[i] = s;
        }
    }

    /**
     * Implicit conversion from vector of convertible element type
     */
    template <typename U>
    HALMD_GPU_ENABLED fixed_vector(fixed_vector<U, 2> const& v,
      typename enable_if<is_convertible<U, dsfloat> >::type* dummy = 0)
    {
        for (size_t i = 0; i < static_size; ++i) {
            (*this)[i] = static_cast<value_type>(v[i]);
        }
    }

    HALMD_GPU_ENABLED fixed_vector(fixed_vector<float, 2> const& v, fixed_vector<float, 2> const& w)
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
      typename enable_if<is_convertible<U, double> >::type* dummy = 0)
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

/**
 * Three-dimensional double precision doubleing-point vector
 */
template <>
struct fixed_vector<double, 3>
  : fixed_array<double, 3>
{
    typedef fixed_array<double, 3> _Base;
    typedef _Base::value_type value_type;
    enum { static_size = _Base::static_size };

    HALMD_GPU_ENABLED fixed_vector() {}

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
    HALMD_GPU_ENABLED explicit fixed_vector(fixed_vector<U, 3> const& v,
      typename enable_if<is_convertible<U, double> >::type* dummy = 0)
    {
        for (size_t i = 0; i < static_size; ++i) {
            (*this)[i] = static_cast<value_type>(v[i]);
        }
    }

#ifdef HALMD_GPU_DOUBLE_PRECISION

    /**
     * Convert from CUDA vector type
     */
    HALMD_GPU_ENABLED fixed_vector(double3 const& v)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
        (*this)[2] = v.z;
    }

    HALMD_GPU_ENABLED fixed_vector(double4 const& v)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
        (*this)[2] = v.z;
    }

    /**
     * Convert to CUDA vector type
     */
    HALMD_GPU_ENABLED operator double3() const
    {
        double3 v;
        v.x = (*this)[0];
        v.y = (*this)[1];
        v.z = (*this)[2];
        return v;
    }

    HALMD_GPU_ENABLED operator double4() const
    {
        double4 v;
        v.x = (*this)[0];
        v.y = (*this)[1];
        v.z = (*this)[2];
        return v;
    }

#endif /* HALMD_GPU_DOUBLE_PRECISION */
};

/**
 * Two-dimensional double precision doubleing-point vector
 */
template <>
struct fixed_vector<double, 2>
  : fixed_array<double, 2>
{
    typedef fixed_array<double, 2> _Base;
    typedef _Base::value_type value_type;
    enum { static_size = _Base::static_size };

    HALMD_GPU_ENABLED fixed_vector() {}

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
    HALMD_GPU_ENABLED explicit fixed_vector(fixed_vector<U, 2> const& v,
      typename enable_if<is_convertible<U, double> >::type* dummy = 0)
    {
        for (size_t i = 0; i < static_size; ++i) {
            (*this)[i] = static_cast<value_type>(v[i]);
        }
    }

#ifdef HALMD_GPU_DOUBLE_PRECISION

    /**
     * Convert from CUDA vector type
     */
    HALMD_GPU_ENABLED fixed_vector(double2 const& v)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
    }

    HALMD_GPU_ENABLED fixed_vector(double3 const& v)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
    }

    HALMD_GPU_ENABLED fixed_vector(double4 const& v)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
    }

    /**
     * Convert to CUDA vector type
     */
    HALMD_GPU_ENABLED operator double2() const
    {
        double2 v;
        v.x = (*this)[0];
        v.y = (*this)[1];
        return v;
    }

    HALMD_GPU_ENABLED operator double3() const
    {
        double3 v;
        v.x = (*this)[0];
        v.y = (*this)[1];
        return v;
    }

    HALMD_GPU_ENABLED operator double4() const
    {
        double4 v;
        v.x = (*this)[0];
        v.y = (*this)[1];
        return v;
    }

#endif /* HALMD_GPU_DOUBLE_PRECISION */
};

HALMD_GPU_USING(algorithm::gpu::tuple, boost::tuple);
HALMD_GPU_USING(algorithm::gpu::tie, boost::tie);
HALMD_GPU_USING(algorithm::gpu::make_tuple, boost::make_tuple);

/**
 * Returns "high" and "low" single precision vector tuple
 */
template <size_t N>
inline HALMD_GPU_ENABLED
tuple<fixed_vector<float, N>, fixed_vector<float, N> >
split(fixed_vector<dsfloat, N> const& v)
{
    fixed_vector<float, N> hi, lo;
    for (size_t i = 0; i < N; ++i) {
        tie(hi[i], lo[i]) = split(v[i]);
    }
    return make_tuple(hi, lo);
}

/**
 * Assignment by elementwise vector addition
 */
template <typename T, typename U, size_t N>
inline HALMD_GPU_ENABLED
typename enable_if<is_convertible<U, T>, fixed_vector<T, N>&>::type
operator+=(fixed_vector<T, N>& v, fixed_vector<U, N> const& w)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] += w[i];
    }
    return v;
}

/**
 * Assignment by elementwise vector subtraction
 */
template <typename T, typename U, size_t N>
inline HALMD_GPU_ENABLED
typename enable_if<is_convertible<U, T>, fixed_vector<T, N>&>::type
operator-=(fixed_vector<T, N>& v, fixed_vector<U, N> const& w)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] -= w[i];
    }
    return v;
}

/**
 * Assignment by scalar multiplication
 */
template <typename T, typename U, size_t N>
inline HALMD_GPU_ENABLED
typename enable_if<is_convertible<U, T>, fixed_vector<T, N>&>::type
operator*=(fixed_vector<T, N>& v, U s)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] *= s;
    }
    return v;
}

/**
 * Assignment by scalar division
 */
template <typename T, typename U, size_t N>
inline HALMD_GPU_ENABLED
typename enable_if<is_convertible<U, T>, fixed_vector<T, N>&>::type
operator/=(fixed_vector<T, N>& v, U s)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] /= s;
    }
    return v;
}

/**
 * Assignment by scalar modulus
 */
template <typename T, typename U, size_t N>
inline HALMD_GPU_ENABLED
typename enable_if<mpl::and_<is_integral<T>, is_integral<U> >, fixed_vector<T, N>&>::type
operator%=(fixed_vector<T, N>& v, U s)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] %= s;
    }
    return v;
}

/**
 * Elementwise vector addition
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<T, N> operator+(fixed_vector<T, N> v, fixed_vector<T, N> const& w)
{
    v += w;
    return v;
}

/**
 * Elementwise vector subtraction
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<T, N> operator-(fixed_vector<T, N> v, fixed_vector<T, N> const& w)
{
    v -= w;
    return v;
}

/**
 * Elementwise change of sign
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<T, N> operator-(fixed_vector<T, N> v)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] = -v[i];
    }
    return v;
}

/**
 * Scalar multiplication
 */
template <typename T, typename U, size_t N>
inline HALMD_GPU_ENABLED
typename enable_if<is_convertible<U, T>, fixed_vector<T, N> >::type
operator*(fixed_vector<T, N> v, U s)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] *= s;
    }
    return v;
}

/**
 * Scalar multiplication
 */
template <typename T, typename U, size_t N>
inline HALMD_GPU_ENABLED
typename enable_if<is_convertible<U, T>, fixed_vector<T, N> >::type
operator*(U s, fixed_vector<T, N> v)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] *= s;
    }
    return v;
}

/**
 * Scalar division
 */
template <typename T, typename U, size_t N>
inline HALMD_GPU_ENABLED
typename enable_if<is_convertible<U, T>, fixed_vector<T, N> >::type
operator/(fixed_vector<T, N> v, U s)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] /= s;
    }
    return v;
}

/**
 * Scalar modulus
 */
template <typename T, typename U, size_t N>
inline HALMD_GPU_ENABLED
typename enable_if<mpl::and_<is_integral<T>, is_integral<U> >, fixed_vector<T, N> >::type
operator%(fixed_vector<T, N> v, U s)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] %= s;
    }
    return v;
}

/**
 * Inner product
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED
T inner_prod(fixed_vector<T, N> const& v, fixed_vector<T, N> const& w)
{
    T s = v[0] * w[0];
    for (size_t i = 1; i < N; ++i) {
        s += v[i] * w[i];
    }
    return s;
}

/**
 * Elementwise vector multiplication
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<T, N> element_prod(fixed_vector<T, N> v, fixed_vector<T, N> const& w)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] *= w[i];
    }
    return v;
}

/**
 * Elementwise vector division
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<T, N> element_div(fixed_vector<T, N> v, fixed_vector<T, N> const& w)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] /= w[i];
    }
    return v;
}

/**
 * Elementwise vector modulus
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED
typename enable_if<is_integral<T>, fixed_vector<T, N> >::type
element_mod(fixed_vector<T, N> v, fixed_vector<T, N> const& w)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] %= w[i];
    }
    return v;
}

/**
 * Elementwise round to nearest integer not greater than argument
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED
typename enable_if<is_floating_point<T>, fixed_vector<T, N> >::type
floor(fixed_vector<T, N> v)
{
    HALMD_GPU_USING(::floor, std::floor);
    for (size_t i = 0; i < N; ++i) {
        v[i] = floor(v[i]);
    }
    return v;
}

/**
 * Elementwise round to nearest integer not less argument
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED
typename enable_if<is_floating_point<T>, fixed_vector<T, N> >::type
ceil(fixed_vector<T, N> v)
{
    HALMD_GPU_USING(::ceil, std::ceil);
    for (size_t i = 0; i < N; ++i) {
        v[i] = ceil(v[i]);
    }
    return v;
}

/**
 * Elementwise square root function
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED
typename enable_if<is_floating_point<T>, fixed_vector<T, N> >::type
sqrt(fixed_vector<T, N> v)
{
    HALMD_GPU_USING(::sqrt, std::sqrt);
    for (size_t i = 0; i < N; ++i) {
        v[i] = sqrt(v[i]);
    }
    return v;
}

/**
 * Elementwise cosine function
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED
typename enable_if<is_floating_point<T>, fixed_vector<T, N> >::type
cos(fixed_vector<T, N> v)
{
    HALMD_GPU_USING(::cos, std::cos);
    for (size_t i = 0; i < N; ++i) {
        v[i] = cos(v[i]);
    }
    return v;
}

/**
 * Elementwise sine function
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED
typename enable_if<is_floating_point<T>, fixed_vector<T, N> >::type
sin(fixed_vector<T, N> v)
{
    HALMD_GPU_USING(::sin, std::sin);
    for (size_t i = 0; i < N; ++i) {
        v[i] = sin(v[i]);
    }
    return v;
}

/**
 * Elementwise absolute value
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED
typename enable_if<is_floating_point<T>, fixed_vector<T, N> >::type
fabs(fixed_vector<T, N> v)
{
    HALMD_GPU_USING(::fabs, std::fabs);
    for (size_t i = 0; i < N; ++i) {
        v[i] = fabs(v[i]);
    }
    return v;
}

/**
 * Floating-point remainder function, round towards nearest integer
 */
template <size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<double, N>
remainder(fixed_vector<double, N> v, fixed_vector<double, N> const& w)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] = ::remainder(v[i], w[i]);
    }
    return v;
}

template <size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<float, N>
remainder(fixed_vector<float, N> v, fixed_vector<float, N> const& w)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] = ::remainderf(v[i], w[i]);
    }
    return v;
}

/**
 * Floating-point remainder function, round towards zero
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED
typename enable_if<is_floating_point<T>, fixed_vector<T, N> >::type
fmod(fixed_vector<T, N> v, fixed_vector<T, N> const& w)
{
    HALMD_GPU_USING(::fmod, std::fmod);
    for (size_t i = 0; i < N; ++i) {
        v[i] = fmod(v[i], w[i]);
    }
    return v;
}

/**
 * Elementwise round to nearest integer
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED
typename enable_if<is_same<T, double>, fixed_vector<T, N> >::type
rint(fixed_vector<T, N> v)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] = ::rint(v[i]);
    }
    return v;
}

template <typename T, size_t N>
inline HALMD_GPU_ENABLED
typename enable_if<is_same<T, float>, fixed_vector<T, N> >::type
rint(fixed_vector<T, N> v)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] = ::rintf(v[i]);
    }
    return v;
}

/**
 * Elementwise round to nearest integer, away from zero
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED
typename enable_if<is_same<T, double>, fixed_vector<T, N> >::type
round(fixed_vector<T, N> v)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] = ::round(v[i]);
    }
    return v;
}

template <typename T, size_t N>
inline HALMD_GPU_ENABLED
typename enable_if<is_same<T, float>, fixed_vector<T, N> >::type
round(fixed_vector<T, N> v)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] = ::roundf(v[i]);
    }
    return v;
}

#ifdef __CUDACC__

/**
 * Fast, accurate floating-point division by s < 2^126
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED
typename enable_if<is_same<T, double>, fixed_vector<T, N> >::type
fdivide(fixed_vector<T, N> v, fixed_vector<T, N> const& w)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] = ::fdivide(v[i], w[i]);
    }
    return v;
}

template <typename T, size_t N>
inline HALMD_GPU_ENABLED
typename enable_if<is_same<T, float>, fixed_vector<T, N> >::type
fdivide(fixed_vector<T, N> v, fixed_vector<T, N> const& w)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] = ::fdividef(v[i], w[i]);
    }
    return v;
}

/**
 * Limit floating-point elements to unit interval [0, 1]
 */
template <typename T, size_t N>
inline HALMD_GPU_ENABLED
typename enable_if<is_same<T, float>, fixed_vector<T, N> >::type
saturate(fixed_vector<T, N> v)
{
    for (size_t i = 0; i < N; ++i) {
        v[i] = ::saturate(v[i]);
    }
    return v;
}

/**
 * Convert floating-point elements to integers, rounding to negative infinity
 */
#ifdef HALMD_GPU_DOUBLE_PRECISION
template <size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<int, N> __double2int_rd(fixed_vector<double, N> const& v)
{
    fixed_vector<int, N> w;
    for (size_t i = 0; i < N; ++i) {
        w[i] = ::__double2int_rd(v[i]);
    }
    return w;
}
#endif

template <size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<int, N> __float2int_rd(fixed_vector<float, N> const& v)
{
    fixed_vector<int, N> w;
    for (size_t i = 0; i < N; ++i) {
        w[i] = ::__float2int_rd(v[i]);
    }
    return w;
}

/**
 * Convert floating-point elements to integers, rounding to nearest even integer
 */
#ifdef HALMD_GPU_DOUBLE_PRECISION
template <size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<int, N> __double2int_rn(fixed_vector<double, N> const& v)
{
    fixed_vector<int, N> w;
    for (size_t i = 0; i < N; ++i) {
        w[i] = ::__double2int_rn(v[i]);
    }
    return w;
}
#endif

template <size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<int, N> __float2int_rn(fixed_vector<float, N> const& v)
{
    fixed_vector<int, N> w;
    for (size_t i = 0; i < N; ++i) {
        w[i] = ::__float2int_rn(v[i]);
    }
    return w;
}

/**
 * Convert floating-point elements to integers, rounding to positive infinity
 */
#ifdef HALMD_GPU_DOUBLE_PRECISION
template <size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<int, N> __double2int_ru(fixed_vector<double, N> const& v)
{
    fixed_vector<int, N> w;
    for (size_t i = 0; i < N; ++i) {
        w[i] = ::__double2int_ru(v[i]);
    }
    return w;
}
#endif

template <size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<int, N> __float2int_ru(fixed_vector<float, N> const& v)
{
    fixed_vector<int, N> w;
    for (size_t i = 0; i < N; ++i) {
        w[i] = ::__float2int_ru(v[i]);
    }
    return w;
}

/**
 * Convert floating-point elements to integers, rounding towards zero
 */
#ifdef HALMD_GPU_DOUBLE_PRECISION
template <size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<int, N> __double2int_rz(fixed_vector<double, N> const& v)
{
    fixed_vector<int, N> w;
    for (size_t i = 0; i < N; ++i) {
        w[i] = ::__double2int_rz(v[i]);
    }
    return w;
}
#endif

template <size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<int, N> __float2int_rz(fixed_vector<float, N> const& v)
{
    fixed_vector<int, N> w;
    for (size_t i = 0; i < N; ++i) {
        w[i] = ::__float2int_rz(v[i]);
    }
    return w;
}

/**
 * Convert floating-point elements to unsigned integers, rounding to negative infinity
 */
#ifdef HALMD_GPU_DOUBLE_PRECISION
template <size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<unsigned int, N> __double2uint_rd(fixed_vector<double, N> const& v)
{
    for (size_t i = 0; i < N; ++i) {
        fixed_vector<unsigned int, N> w;
        w[i] = ::__double2uint_rd(v[i]);
    }
    return w;
}
#endif

template <size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<unsigned int, N> __float2uint_rd(fixed_vector<float, N> const& v)
{
    fixed_vector<unsigned int, N> w;
    for (size_t i = 0; i < N; ++i) {
        w[i] = ::__float2uint_rd(v[i]);
    }
    return w;
}

/**
 * Convert floating-point elements to unsigned integers, rounding to nearest even integer
 */
#ifdef HALMD_GPU_DOUBLE_PRECISION
template <size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<unsigned int, N> __double2uint_rn(fixed_vector<double, N> const& v)
{
    fixed_vector<unsigned int, N> w;
    for (size_t i = 0; i < N; ++i) {
        w[i] = ::__double2uint_rn(v[i]);
    }
    return w;
}
#endif

template <size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<unsigned int, N> __float2uint_rn(fixed_vector<float, N> const& v)
{
    fixed_vector<unsigned int, N> w;
    for (size_t i = 0; i < N; ++i) {
        w[i] = ::__float2uint_rn(v[i]);
    }
    return w;
}

/**
 * Convert floating-point elements to unsigned integers, rounding to positive infinity
 */
#ifdef HALMD_GPU_DOUBLE_PRECISION
template <size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<unsigned int, N> __double2uint_ru(fixed_vector<double, N> const& v)
{
    fixed_vector<unsigned int, N> w;
    for (size_t i = 0; i < N; ++i) {
        w[i] = ::__double2uint_ru(v[i]);
    }
    return w;
}
#endif

template <size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<unsigned int, N> __float2uint_ru(fixed_vector<float, N> const& v)
{
    fixed_vector<unsigned int, N> w;
    for (size_t i = 0; i < N; ++i) {
        w[i] = ::__float2uint_ru(v[i]);
    }
    return w;
}

/**
 * Convert floating-point elements to unsigned integers, rounding towards zero
 */
#ifdef HALMD_GPU_DOUBLE_PRECISION
template <size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<unsigned int, N> __double2uint_rz(fixed_vector<double, N> const& v)
{
    fixed_vector<unsigned int, N> w;
    for (size_t i = 0; i < N; ++i) {
        w[i] = ::__double2uint_rz(v[i]);
    }
    return w;
}
#endif

template <size_t N>
inline HALMD_GPU_ENABLED
fixed_vector<unsigned int, N> __float2uint_rz(fixed_vector<float, N> const& v)
{
    fixed_vector<unsigned int, N> w;
    for (size_t i = 0; i < N; ++i) {
        w[i] = ::__float2uint_rz(v[i]);
    }
    return w;
}

#else /* ! __CUDACC__ */

/**
 * Write vector elements to output stream
 */
template <typename T, size_t N>
inline std::ostream& operator<<(std::ostream& os, fixed_vector<T, N> const& v)
{
    os << v[0];
    for (size_t i = 1; i < N; ++i) {
        os << " " << v[i];
    }
    return os;
}

/**
 * Read vector elements from input stream
 */
template <typename T, size_t N>
inline std::istream& operator>>(std::istream& is, fixed_vector<T, N>& v)
{
    for (size_t i = 0; i < N; ++i) {
        is >> v[i];
    }
    return is;
}

#endif /* ! __CUDACC__ */

}}} // namespace detail::numeric::blas

// import into top-level namespace
using detail::numeric::blas::fixed_vector;

} // namespace halmd

#endif /* ! HALMD_NUMERIC_BLAS_FIXED_VECTOR_HPP */
