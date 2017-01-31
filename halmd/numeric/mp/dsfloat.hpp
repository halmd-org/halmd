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

#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/type_traits/is_floating_point.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>

#include <halmd/numeric/mp/dsfun.hpp>
#include <halmd/utility/tuple.hpp>
#include <cuda_wrapper/cuda_wrapper.hpp>

namespace halmd {
namespace detail {
namespace numeric {
namespace mp {

struct dsfloat;
template<typename T>
struct dsfloat_ptr;
template<typename T>
struct dsfloat_const_ptr;
#ifndef __CUDACC__
template<typename T>
class dsfloat_vector;
#endif

} // namespace mp
} // namespace numeric
} // namespace detail

// import into top-level namespace
using detail::numeric::mp::dsfloat;
using detail::numeric::mp::dsfloat_ptr;
using detail::numeric::mp::dsfloat_const_ptr;
#ifndef __CUDACC__
using detail::numeric::mp::dsfloat_vector;
#endif

} // namespace halmd

namespace boost
{

//
// In order to make dsfloat completely act like a number, it is necessary
// to specify it as floating point number by the means of
// boost::is_floating_point<>. Therefore, add a specialization for
// boost::is_floating_point<dsfloat>.
//
template<>
struct is_floating_point<halmd::dsfloat> : public boost::true_type {};

} // namespace boost

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
struct dsfloat_ptr {
#ifndef __CUDACC__
    explicit dsfloat_ptr(dsfloat_vector<T>& v);
#endif

    T* hi;
    T* lo;

    HALMD_GPU_ENABLED halmd::tuple<T&, T&> operator[] (unsigned int idx) {
        return tie(hi[idx], lo[idx]);
    };
    HALMD_GPU_ENABLED halmd::tuple<T const&, T const&> operator[] (unsigned int idx) const {
        return tie(hi[idx], lo[idx]);
    };
};

template<typename T>
struct dsfloat_const_ptr {
#ifndef __CUDACC__
    explicit dsfloat_const_ptr(dsfloat_vector<T> const& v);
#endif

    T const* hi;
    T const* lo;

    HALMD_GPU_ENABLED halmd::tuple<T const&, T const&> operator[] (unsigned int idx) const {
        return tie(hi[idx], lo[idx]);
    };
};

#ifndef __CUDACC__
template<typename T>
class dsfloat_vector {
public:
    typedef dsfloat_vector<T> vector_type;
    typedef T value_type;
    typedef dsfloat_ptr<T> pointer;
    typedef dsfloat_const_ptr<T> const const_pointer;
    typedef size_t size_type;

    dsfloat_vector(size_type size) : data_(size * 2)
    {
    }

    size_type size() const
    {
        return data_.size() / 2;
    }

    size_type capacity() const
    {
        return data_.capacity() / 2;
    }

    void resize(size_type size)
    {
        data_.resize(size * 2);
    }

    void reserve(size_type size)
    {
        data_.reserve(size * 2);
    }

    void swap(vector_type& v)
    {
        data_.swap(v.data_);
    }

    pointer data()
    {
        return pointer(*this);
    }

    operator pointer()
    {
        return data();
    }

    const_pointer data() const
    {
        return const_pointer(*this);
    }

    operator const_pointer() const
    {
        return data();
    }

    cuda::vector<T>& storage() {
        return data_;
    }

    cuda::vector<T> const& storage() const {
        return data_;
    }

private:
    cuda::vector<T> data_;
};

template<typename T>
inline dsfloat_ptr<T>::dsfloat_ptr(dsfloat_vector<T>& v)
        : hi(&*(v.storage().begin()))
        , lo(&*(v.storage().begin()+v.size()))
{
}

template<typename T>
inline dsfloat_const_ptr<T>::dsfloat_const_ptr(dsfloat_vector<T> const& v)
  : hi(&*(v.storage().begin()))
    , lo(&*(v.storage().begin()+v.size()))
{
}
#endif

} // namespace mp
} // namespace numeric
} // namespace detail
} // namespace halmd

#endif /* ! HALMD_NUMERIC_MP_DSFLOAT_CUH */
