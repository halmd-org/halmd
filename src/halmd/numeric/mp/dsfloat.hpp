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

#ifndef HALMD_NUMERIC_MP_DSFLOAT_CUH
#define HALMD_NUMERIC_MP_DSFLOAT_CUH

#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>

#include <halmd/numeric/mp/dsfun.hpp>
#ifdef __CUDACC__
# include <halmd/algorithm/gpu/tuple.cuh>
#else
# include <boost/tuple/tuple.hpp>
#endif

namespace halmd
{
namespace detail { namespace numeric { namespace mp
{

#ifdef __CUDACC__
using algorithm::gpu::tuple;
using algorithm::gpu::make_tuple;
using algorithm::gpu::tie;
#else
using boost::tuple;
using boost::make_tuple;
using boost::tie;
#endif

/**
 * Double-single floating point value
 */
struct dsfloat
{
    float hi, lo;

    HALMD_GPU_ENABLED dsfloat()
    {}

    HALMD_GPU_ENABLED dsfloat(float a0, float a1)
      : hi(a0), lo(a1)
    {}

    template <typename T>
    HALMD_GPU_ENABLED dsfloat(T a0,
      typename boost::enable_if<boost::is_same<T, float> >::type* dummy = 0)
    {
        dsfeq(hi, lo, a0);
    }

    template <typename T>
    HALMD_GPU_ENABLED dsfloat(T a0,
      typename boost::enable_if<boost::is_integral<T> >::type* dummy = 0)
    {
        dsfeq(hi, lo, a0);
    }

    template <typename T>
    HALMD_GPU_ENABLED dsfloat(T a,
      typename boost::enable_if<boost::is_same<T, double> >::type* dummy = 0)
    {
        dsdeq(hi, lo, a);
    }

    /**
     * Returns "high" single precision floating-point value
     */
    HALMD_GPU_ENABLED operator float() const
    {
        return hi;
    }

    /**
     * Returns double precision value if supported natively
     */
    HALMD_GPU_ENABLED operator double() const
    {
        return static_cast<double>(hi) + lo;
    }
};

/**
 * Returns "high" and "low" single precision floating-point tuple
 */
inline HALMD_GPU_ENABLED tuple<float, float> split(dsfloat const& v)
{
    return make_tuple(v.hi, v.lo);
}

/**
 * Addition by assignment
 */
inline HALMD_GPU_ENABLED dsfloat& operator+=(dsfloat& v, dsfloat const& w)
{
    dsadd(v.hi, v.lo, v.hi, v.lo, w.hi, w.lo);
    return v;
}

/**
 * Subtraction by assignment
 */
inline HALMD_GPU_ENABLED dsfloat& operator-=(dsfloat& v, dsfloat const& w)
{
    dssub(v.hi, v.lo, v.hi, v.lo, w.hi, w.lo);
    return v;
}

/**
 * Multiplication by assignment
 */
inline HALMD_GPU_ENABLED dsfloat& operator*=(dsfloat& v, dsfloat const& w)
{
    dsmul(v.hi, v.lo, v.hi, v.lo, w.hi, w.lo);
    return v;
}

/**
 * Division by assignment
 */
inline HALMD_GPU_ENABLED dsfloat& operator/=(dsfloat& v, dsfloat const& w)
{
    dsdiv(v.hi, v.lo, v.hi, v.lo, w.hi, w.lo);
    return v;
}

/**
 * Addition
 */
inline HALMD_GPU_ENABLED dsfloat operator+(dsfloat v, dsfloat const& w)
{
    v += w;
    return v;
}

template <typename T>
inline HALMD_GPU_ENABLED typename boost::enable_if<boost::is_arithmetic<T>, dsfloat>::type
operator+(T const& v, dsfloat const& w)
{
    return static_cast<dsfloat>(v) + w;
}

template <typename T>
inline HALMD_GPU_ENABLED typename boost::enable_if<boost::is_arithmetic<T>, dsfloat>::type
operator+(dsfloat const& v, T const& w)
{
    return v + static_cast<dsfloat>(w);
}

/**
 * Subtraction
 */
inline HALMD_GPU_ENABLED dsfloat operator-(dsfloat v, dsfloat const& w)
{
    v -= w;
    return v;
}

template <typename T>
inline HALMD_GPU_ENABLED typename boost::enable_if<boost::is_arithmetic<T>, dsfloat>::type
operator-(T const& v, dsfloat const& w)
{
    return static_cast<dsfloat>(v) - w;
}

template <typename T>
inline HALMD_GPU_ENABLED typename boost::enable_if<boost::is_arithmetic<T>, dsfloat>::type
operator-(dsfloat const& v, T const& w)
{
    return v - static_cast<dsfloat>(w);
}

/**
 * Multiplication
 */
inline HALMD_GPU_ENABLED dsfloat operator*(dsfloat v, dsfloat const& w)
{
    v *= w;
    return v;
}

template <typename T>
inline HALMD_GPU_ENABLED typename boost::enable_if<boost::is_arithmetic<T>, dsfloat>::type
operator*(T const& v, dsfloat const& w)
{
    return static_cast<dsfloat>(v) * w;
}

template <typename T>
inline HALMD_GPU_ENABLED typename boost::enable_if<boost::is_arithmetic<T>, dsfloat>::type
operator*(dsfloat const& v, T const& w)
{
    return v * static_cast<dsfloat>(w);
}

inline HALMD_GPU_ENABLED dsfloat operator/(dsfloat v, dsfloat const& w)
{
    v /= w;
    return v;
}

/**
 * Division
 */
template <typename T>
inline HALMD_GPU_ENABLED typename boost::enable_if<boost::is_arithmetic<T>, dsfloat>::type
operator/(T const& v, dsfloat const& w)
{
    return static_cast<dsfloat>(v) / w;
}

template <typename T>
inline HALMD_GPU_ENABLED typename boost::enable_if<boost::is_arithmetic<T>, dsfloat>::type
operator/(dsfloat const& v, T const& w)
{
    return v / static_cast<dsfloat>(w);
}

/**
 * Square root function
 */
inline HALMD_GPU_ENABLED dsfloat sqrt(dsfloat v)
{
    dsfloat w;
    dssqrt(w.hi, w.lo, v.hi, v.lo);
    return w;
}

}}} // namespace detail::numeric::mp

// import into top-level namespace
using detail::numeric::mp::dsfloat;

} // namespace halmd

#endif /* ! HALMD_NUMERIC_MP_DSFLOAT_CUH */
