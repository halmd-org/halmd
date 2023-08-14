/*
 * Copyright © 2008-2013 Felix Höfling
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

#ifndef HALMD_NUMERIC_BLAS_DETAIL_SIZE_6_HPP
#define HALMD_NUMERIC_BLAS_DETAIL_SIZE_6_HPP

#include <halmd/config.hpp>

#include <algorithm>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/utility/enable_if.hpp>
#include <cassert>
#include <initializer_list>
#include <type_traits>

// CUDA vector types for host compiler
#ifdef HALMD_WITH_GPU
/* Disable warning for CUDA 5.5 headers emitted by Clang:
 *   /usr/local/cuda-5.5/include/cuda_runtime.h:225:33: warning: function
 *   'cudaMallocHost' is not needed and will not be emitted [-Wunneeded-internal-declaration]
 */
# if (defined(__clang__) && __clang_major__ >= 3 && __clang_minor__ > 2)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wunneeded-internal-declaration"
#  include <cuda_runtime.h>
#  pragma GCC diagnostic pop
# else
#  include <cuda_runtime.h>
# endif
#endif

#include <halmd/numeric/blas/detail/array.hpp>
#include <halmd/numeric/blas/detail/vector.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>
#include <halmd/utility/tuple.hpp>

namespace halmd {
namespace numeric {
namespace blas {
namespace detail {

/**
 * Six-dimensional single precision floating-point vector
 */
template <>
struct fixed_vector<float, 6>
  : fixed_array<float, 6>
{
    typedef fixed_array<float, 6> _Base;
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
    HALMD_GPU_ENABLED explicit fixed_vector(fixed_vector<U, 6> const& v,
      typename boost::enable_if<boost::is_convertible<U, float> >::type* dummy = 0)
    {
        for (size_t i = 0; i < static_size; ++i) {
            (*this)[i] = static_cast<value_type>(v[i]);
        }
    }

#ifdef HALMD_WITH_GPU

    /**
     * Convert from CUDA vector type
     */
    HALMD_GPU_ENABLED fixed_vector(float4 const& v, float2 const& w)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
        (*this)[2] = v.z;
        (*this)[3] = v.w;
        (*this)[4] = w.x;
        (*this)[5] = w.y;
    }

    HALMD_GPU_ENABLED fixed_vector(float3 const& v, float3 const& w)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
        (*this)[2] = v.z;
        (*this)[3] = w.x;
        (*this)[4] = w.y;
        (*this)[5] = w.z;
    }

    /**
     * Convert to CUDA vector types
     */
    HALMD_GPU_ENABLED operator tuple<float4, float2>() const
    {
        float4 v;
        float2 w;
        v.x = (*this)[0];
        v.y = (*this)[1];
        v.z = (*this)[2];
        v.w = (*this)[3];
        w.x = (*this)[4];
        w.y = (*this)[5];
        return make_tuple(v, w);
    }

    HALMD_GPU_ENABLED operator tuple<float3, float3>() const
    {
        float3 v;
        float3 w;
        v.x = (*this)[0];
        v.y = (*this)[1];
        v.z = (*this)[2];
        w.x = (*this)[3];
        w.y = (*this)[4];
        w.z = (*this)[5];
        return make_tuple(v, w);
    }

#endif /* HALMD_WITH_GPU */
};

/**
 * Six-dimensional unsigned integer vector
 */
template <>
struct fixed_vector<unsigned int, 6>
  : fixed_array<unsigned int, 6>
{
    typedef fixed_array<unsigned int, 6> _Base;
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
    HALMD_GPU_ENABLED explicit fixed_vector(fixed_vector<U, 6> const& v,
      typename boost::enable_if<boost::is_convertible<U, unsigned int> >::type* dummy = 0)
    {
        for (size_t i = 0; i < static_size; ++i) {
            (*this)[i] = static_cast<value_type>(v[i]);
        }
    }

#ifdef HALMD_WITH_GPU

    /**
     * Convert from CUDA vector type
     */
    HALMD_GPU_ENABLED fixed_vector(uint4 const& v, uint2 const& w)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
        (*this)[2] = v.z;
        (*this)[3] = v.w;
        (*this)[4] = w.x;
        (*this)[5] = w.y;
    }

    HALMD_GPU_ENABLED fixed_vector(uint3 const& v, uint3 const& w)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
        (*this)[2] = v.z;
        (*this)[3] = w.x;
        (*this)[4] = w.y;
        (*this)[5] = w.z;
    }

    /**
     * Convert to CUDA vector types
     */
    HALMD_GPU_ENABLED operator tuple<uint4, uint2>() const
    {
        uint4 v;
        uint2 w;
        v.x = (*this)[0];
        v.y = (*this)[1];
        v.z = (*this)[2];
        v.w = (*this)[3];
        w.x = (*this)[4];
        w.y = (*this)[5];
        return make_tuple(v, w);
    }

    HALMD_GPU_ENABLED operator tuple<uint3, uint3>() const
    {
        uint3 v;
        uint3 w;
        v.x = (*this)[0];
        v.y = (*this)[1];
        v.z = (*this)[2];
        w.x = (*this)[3];
        w.y = (*this)[4];
        w.z = (*this)[5];
        return make_tuple(v, w);
    }

#endif /* HALMD_WITH_GPU */
};

/**
 * Six-dimensional integer vector
 */
template <>
struct fixed_vector<int, 6>
  : fixed_array<int, 6>
{
    typedef fixed_array<int, 6> _Base;
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
    HALMD_GPU_ENABLED explicit fixed_vector(fixed_vector<U, 6> const& v,
      typename boost::enable_if<boost::is_convertible<U, int> >::type* dummy = 0)
    {
        for (size_t i = 0; i < static_size; ++i) {
            (*this)[i] = static_cast<value_type>(v[i]);
        }
    }

#ifdef HALMD_WITH_GPU

    /**
     * Convert from CUDA vector type
     */
    HALMD_GPU_ENABLED fixed_vector(int4 const& v, int2 const& w)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
        (*this)[2] = v.z;
        (*this)[3] = v.w;
        (*this)[4] = w.x;
        (*this)[5] = w.y;
    }

    HALMD_GPU_ENABLED fixed_vector(int3 const& v, int3 const& w)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
        (*this)[2] = v.z;
        (*this)[3] = w.x;
        (*this)[4] = w.y;
        (*this)[5] = w.z;
    }

    /**
     * Convert to CUDA vector types
     */
    HALMD_GPU_ENABLED operator tuple<int4, int2>() const
    {
        int4 v;
        int2 w;
        v.x = (*this)[0];
        v.y = (*this)[1];
        v.z = (*this)[2];
        v.w = (*this)[3];
        w.x = (*this)[4];
        w.y = (*this)[5];
        return make_tuple(v, w);
    }

    HALMD_GPU_ENABLED operator tuple<int3, int3>() const
    {
        int3 v;
        int3 w;
        v.x = (*this)[0];
        v.y = (*this)[1];
        v.z = (*this)[2];
        w.x = (*this)[3];
        w.y = (*this)[4];
        w.z = (*this)[5];
        return make_tuple(v, w);
    }

#endif /* HALMD_WITH_GPU */
};

/**
 * Six-dimensional double-single precision floating-point vector
 */
template <>
struct fixed_vector<dsfloat, 6> : fixed_array<dsfloat, 6>
{
    typedef fixed_array<dsfloat, 6> _Base;
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
      typename boost::enable_if<boost::is_convertible<U, dsfloat> >::type* dummy = 0)
    {
        for (size_t i = 0; i < static_size; ++i) {
            (*this)[i] = s;
        }
    }

    /**
     * Implicit conversion from vector of convertible element type
     */
    template <typename U>
    HALMD_GPU_ENABLED fixed_vector(fixed_vector<U, 6> const& v,
      typename boost::enable_if<boost::is_convertible<U, dsfloat> >::type* dummy = 0)
    {
        for (size_t i = 0; i < static_size; ++i) {
            (*this)[i] = static_cast<value_type>(v[i]);
        }
    }

    HALMD_GPU_ENABLED fixed_vector(fixed_vector<float, 6> const& v, fixed_vector<float, 6> const& w)
    {
        for (size_t i = 0; i < static_size; ++i) {
            (*this)[i] = value_type(v[i], w[i]);
        }
    }
};

/**
 * Six-dimensional double precision floating-point vector
 */
template <>
struct fixed_vector<double, 6>
  : fixed_array<double, 6>
{
    typedef fixed_array<double, 6> _Base;
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
    HALMD_GPU_ENABLED explicit fixed_vector(fixed_vector<U, 6> const& v,
      typename boost::enable_if<boost::is_convertible<U, double> >::type* dummy = 0)
    {
        for (size_t i = 0; i < static_size; ++i) {
            (*this)[i] = static_cast<value_type>(v[i]);
        }
    }

#ifdef HALMD_GPU_DOUBLE_PRECISION

    /**
     * Convert from CUDA vector type
     */
    HALMD_GPU_ENABLED fixed_vector(double4 const& v, double2 const& w)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
        (*this)[2] = v.z;
        (*this)[3] = v.w;
        (*this)[4] = w.x;
        (*this)[5] = w.y;
    }

    HALMD_GPU_ENABLED fixed_vector(double3 const& v, double3 const& w)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
        (*this)[2] = v.z;
        (*this)[3] = w.x;
        (*this)[4] = w.y;
        (*this)[5] = w.z;
    }

    /**
     * Convert to CUDA vector types
     */
    HALMD_GPU_ENABLED operator tuple<double4, double2>() const
    {
        double4 v;
        double2 w;
        v.x = (*this)[0];
        v.y = (*this)[1];
        v.z = (*this)[2];
        v.w = (*this)[3];
        w.x = (*this)[4];
        w.y = (*this)[5];
        return make_tuple(v, w);
    }

    HALMD_GPU_ENABLED operator tuple<double3, double3>() const
    {
        double3 v;
        double3 w;
        v.x = (*this)[0];
        v.y = (*this)[1];
        v.z = (*this)[2];
        w.x = (*this)[3];
        w.y = (*this)[4];
        w.z = (*this)[5];
        return make_tuple(v, w);
    }

#endif /* HALMD_GPU_DOUBLE_PRECISION */
};

#ifdef __CUDACC__

/**
 * Split in tuple of coalesced types
 */
inline HALMD_GPU_ENABLED
tuple<float4, float2> split(fixed_vector<float, 6> const& v)
{
    return static_cast<tuple<float4, float2> >(v);
}

inline HALMD_GPU_ENABLED
tuple<uint4, uint2> split(fixed_vector<uint, 6> const& v)
{
    return static_cast<tuple<uint4, uint2> >(v);
}

inline HALMD_GPU_ENABLED
tuple<int4, int2> split(fixed_vector<int, 6> const& v)
{
    return static_cast<tuple<int4, int2> >(v);
}

#ifdef HALMD_GPU_DOUBLE_PRECISION

inline HALMD_GPU_ENABLED
tuple<double4, double2> split(fixed_vector<double, 6> const& v)
{
    return static_cast<tuple<double4, double2> >(v);
}

#endif /* HALMD_GPU_DOUBLE_PRECISION */

#endif  // __CUDACC__

} // namespace detail
} // namespace blas
} // namespace numeric
} // namespace halmd

#endif /* ! HALMD_NUMERIC_BLAS_DETAIL_SIZE_6_HPP */
