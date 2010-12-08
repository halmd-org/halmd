/*
 * Copyright © 2008-2010  Peter Colberg and Felix Höfling
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

#ifndef HALMD_NUMERIC_BLAS_FIXED_VECTOR_SIZE_3_HPP
#define HALMD_NUMERIC_BLAS_FIXED_VECTOR_SIZE_3_HPP

#include <boost/type_traits/is_convertible.hpp>
#include <boost/utility/enable_if.hpp>
#ifdef WITH_CUDA
# include <cuda_runtime.h> // CUDA vector types for host compiler
#endif

#include <halmd/config.hpp>
#include <halmd/numeric/blas/fixed_array.hpp>
#include <halmd/numeric/blas/fixed_vector/size_N.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>
#ifdef __CUDACC__
# include <halmd/algorithm/gpu/tuple.cuh>
#else
# include <boost/tuple/tuple.hpp>
#endif

namespace halmd
{
namespace detail { namespace numeric { namespace blas
{

#ifdef WITH_CUDA
HALMD_GPU_USING(algorithm::gpu::tuple, boost::tuple);
HALMD_GPU_USING(algorithm::gpu::make_tuple, boost::make_tuple);
#endif

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
      typename boost::enable_if<boost::is_convertible<U, float> >::type* dummy = 0)
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

    HALMD_GPU_ENABLED fixed_vector(float2 const& v, float1 w)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
        (*this)[2] = w.x;
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

    HALMD_GPU_ENABLED operator tuple<float2, float1>() const
    {
        float2 v;
        float1 w;
        v.x = (*this)[0];
        v.y = (*this)[1];
        w.x = (*this)[2];
        return make_tuple(v, w);
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
      typename boost::enable_if<boost::is_convertible<U, unsigned int> >::type* dummy = 0)
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

    HALMD_GPU_ENABLED fixed_vector(uint2 const& v, uint1 w)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
        (*this)[2] = w.x;
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

    HALMD_GPU_ENABLED operator tuple<uint2, uint1>() const
    {
        uint2 v;
        uint1 w;
        v.x = (*this)[0];
        v.y = (*this)[1];
        w.x = (*this)[2];
        return make_tuple(v, w);
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
      typename boost::enable_if<boost::is_convertible<U, int> >::type* dummy = 0)
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

    HALMD_GPU_ENABLED fixed_vector(int2 const& v, int1 w)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
        (*this)[2] = w.x;
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

    HALMD_GPU_ENABLED operator tuple<int2, int1>() const
    {
        int2 v;
        int1 w;
        v.x = (*this)[0];
        v.y = (*this)[1];
        w.x = (*this)[2];
        return make_tuple(v, w);
    }

#endif /* WITH_CUDA */
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
    HALMD_GPU_ENABLED fixed_vector(fixed_vector<U, 3> const& v,
      typename boost::enable_if<boost::is_convertible<U, dsfloat> >::type* dummy = 0)
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
 * Three-dimensional double precision floating-point vector
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

    HALMD_GPU_ENABLED fixed_vector(double2 const& v, double1 w)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
        (*this)[2] = w.x;
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

    HALMD_GPU_ENABLED operator tuple<double2, double1>() const
    {
        double2 v;
        double1 w;
        v.x = (*this)[0];
        v.y = (*this)[1];
        w.x = (*this)[2];
        return make_tuple(v, w);
    }

#endif /* HALMD_GPU_DOUBLE_PRECISION */
};

#ifdef __CUDACC__

/**
 * Split in tuple of coalesced types
 */
inline HALMD_GPU_ENABLED
tuple<float2, float1> split(fixed_vector<float, 3> const& v)
{
    return static_cast<tuple<float2, float1> >(v);
}

inline HALMD_GPU_ENABLED
tuple<uint2, uint1> split(fixed_vector<unsigned int, 3> const& v)
{
    return static_cast<tuple<uint2, uint1> >(v);
}

inline HALMD_GPU_ENABLED
tuple<int2, int1> split(fixed_vector<int, 3> const& v)
{
    return static_cast<tuple<int2, int1> >(v);
}

#ifdef HALMD_GPU_DOUBLE_PRECISION

inline HALMD_GPU_ENABLED
tuple<double2, double1> split(fixed_vector<double, 3> const& v)
{
    return static_cast<tuple<double2, double1> >(v);
}

#endif /* HALMD_GPU_DOUBLE_PRECISION */

#endif  // __CUDACC__

}}} // namespace detail::numeric::blas

} // namespace halmd

#endif /* ! HALMD_NUMERIC_BLAS_FIXED_VECTOR_SIZE_3_HPP */
