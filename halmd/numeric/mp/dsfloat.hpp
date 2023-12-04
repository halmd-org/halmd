/*
 * Copyright © 2008-2010  Peter Colberg
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

#ifndef HALMD_NUMERIC_MP_DSFLOAT_CUH
#define HALMD_NUMERIC_MP_DSFLOAT_CUH

#include <halmd/numeric/mp/dsfun.hpp>
#include <halmd/utility/tuple.hpp>
#ifdef HALMD_WITH_GPU
# include <cuda_wrapper/cuda_wrapper.hpp>
#endif

#ifndef __CUDACC__
# include <ostream>
#endif
#include <type_traits>

namespace halmd {
namespace detail {
namespace numeric {
namespace mp {

// forward declaration
struct dsfloat;

} // namespace mp
} // namespace numeric
} // namespace detail
} // namespace halmd

namespace std {

//
// In order to make dsfloat completely act like a number, it is necessary
// to specify it as floating point number by the means of
// std::is_floating_point<>. Therefore, add a specialization for
// std::is_floating_point<dsfloat>.
//
template<>
struct is_floating_point<halmd::detail::numeric::mp::dsfloat> : public std::true_type {};

} // namespace std

namespace halmd {
namespace detail {
namespace numeric {
namespace mp {

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
      typename std::enable_if<std::is_same<T, float>::value>::type* dummy = 0)
    {
        dsfeq(hi, lo, a0);
    }

    template <typename T>
    HALMD_GPU_ENABLED dsfloat(T a0,
      typename std::enable_if<std::is_integral<T>::value>::type* dummy = 0)
    {
        dsfeq(hi, lo, a0);
    }

    template <typename T>
    HALMD_GPU_ENABLED dsfloat(T a,
      typename std::enable_if<std::is_same<T, double>::value>::type* dummy = 0)
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

    HALMD_GPU_ENABLED bool operator<(dsfloat const& rhs) const {
        if (hi < rhs.hi) return true;
        if (rhs.hi < hi) return false;
        return lo < rhs.lo;
    }
    HALMD_GPU_ENABLED bool operator>(dsfloat const& rhs) const {
        return rhs.operator<(*this);
    }
    HALMD_GPU_ENABLED bool operator<=(dsfloat const& rhs) const {
        return !rhs.operator<(*this);
    }
    HALMD_GPU_ENABLED bool operator>=(dsfloat const& rhs) const {
        return !operator<(rhs);
    }
    HALMD_GPU_ENABLED bool operator==(dsfloat const& rhs) const {
        return (hi == rhs.hi) && (lo == rhs.lo);
    }
    HALMD_GPU_ENABLED bool operator!=(dsfloat const& rhs) const {
        return !operator==(rhs);
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
 * Sign
 */
inline HALMD_GPU_ENABLED dsfloat operator-(dsfloat v)
{
    v.hi = -v.hi;
    v.lo = -v.lo;
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
inline HALMD_GPU_ENABLED typename std::enable_if<std::is_arithmetic<T>::value, dsfloat>::type
operator+(T const& v, dsfloat const& w)
{
    return static_cast<dsfloat>(v) + w;
}

template <typename T>
inline HALMD_GPU_ENABLED typename std::enable_if<std::is_arithmetic<T>::value, dsfloat>::type
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
inline HALMD_GPU_ENABLED typename std::enable_if<std::is_arithmetic<T>::value, dsfloat>::type
operator-(T const& v, dsfloat const& w)
{
    return static_cast<dsfloat>(v) - w;
}

template <typename T>
inline HALMD_GPU_ENABLED typename std::enable_if<std::is_arithmetic<T>::value, dsfloat>::type
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
inline HALMD_GPU_ENABLED typename std::enable_if<std::is_arithmetic<T>::value, dsfloat>::type
operator*(T const& v, dsfloat const& w)
{
    return static_cast<dsfloat>(v) * w;
}

template <typename T>
inline HALMD_GPU_ENABLED typename std::enable_if<std::is_arithmetic<T>::value, dsfloat>::type
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
inline HALMD_GPU_ENABLED typename std::enable_if<std::is_arithmetic<T>::value, dsfloat>::type
operator/(T const& v, dsfloat const& w)
{
    return static_cast<dsfloat>(v) / w;
}

template <typename T>
inline HALMD_GPU_ENABLED typename std::enable_if<std::is_arithmetic<T>::value, dsfloat>::type
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

/**
 * Maximum value
 */
inline HALMD_GPU_ENABLED dsfloat max(dsfloat const& v, dsfloat const& w)
{
    return v.hi == w.hi ? (v.lo >= w.lo ? v : w) : (v.hi > w.hi ? v : w);
}

/**
 * Minimum value
 */
inline HALMD_GPU_ENABLED dsfloat min(dsfloat const& v, dsfloat const& w)
{
    return v.hi == w.hi ? (v.lo <= w.lo ? v : w) : (v.hi < w.hi ? v : w);
}

template<typename T>
struct dsfloat_ptr
{
    T* hi;
    T* lo;

    HALMD_GPU_ENABLED halmd::tuple<T&, T&> operator[] (unsigned int idx)
    {
        return tie(hi[idx], lo[idx]);
    };
    HALMD_GPU_ENABLED halmd::tuple<T const&, T const&> operator[] (unsigned int idx) const
    {
        return tie(hi[idx], lo[idx]);
    };

    operator T*() const
    {
        return hi;
    }
};

template<typename T>
struct dsfloat_const_ptr
{
    T const* hi;
    T const* lo;

    HALMD_GPU_ENABLED halmd::tuple<T const&, T const&> operator[] (unsigned int idx) const
    {
        return tie(hi[idx], lo[idx]);
    };

    operator T const*() const
    {
        return hi;
    }
};

#ifndef __CUDACC__
inline std::ostream& operator<<(std::ostream& p, dsfloat const& val) {
    p << double(val);
    return p;
}
#endif

} // namespace mp
} // namespace numeric
} // namespace detail

// import into top-level namespace
using detail::numeric::mp::dsfloat;
using detail::numeric::mp::dsfloat_ptr;
using detail::numeric::mp::dsfloat_const_ptr;

} // namespace halmd

#endif /* ! HALMD_NUMERIC_MP_DSFLOAT_CUH */
