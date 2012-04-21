/*
 * Copyright © 2008-2012  Peter Colberg
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

#ifndef HALMD_NUMERIC_BLAS_FIXED_VECTOR_CUDA_VECTOR_CONVERTER_HPP
#define HALMD_NUMERIC_BLAS_FIXED_VECTOR_CUDA_VECTOR_CONVERTER_HPP

#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>

#include <halmd/config.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>
#include <halmd/utility/tuple.hpp>

namespace halmd {
namespace detail {
namespace numeric {
namespace blas {

#ifndef __CUDACC__

/**
 * Returns bit representation of an int value stored in a float value.
 *
 * This function uses a union cast to convert between float and int.
 */
inline float __int_as_float(int value)
{
    union { int i; float f; } volatile u;
    u.i = value;
    return u.f;
}

/**
 * Returns bit representation of a float value stored in an int value.
 *
 * This function uses a union cast to convert between float and int.
 */
inline int __float_as_int(float value)
{
    union { float f; int i; } volatile u;
    u.f = value;
    return u.i;
}

#endif /* ! __CUDACC__ */

template <typename L, typename R>
struct cuda_vector_converter;

/**
 * Pack tuple of single-precision floating-point value and scalar into CUDA vector.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::is_same<T, typename cuda_vector_converter<T, tuple<U const, V const> >::result_type>, void>::type
operator<<=(T& left, tuple<U&, V&> right)
{
    left = cuda_vector_converter<T, tuple<U const, V const> >::apply(get<0>(right), get<1>(right));
}

/**
 * Pack tuple of single-precision floating-point value and scalar into CUDA vector.
 *
 * This function accepts a tuple of two CUDA vectors for signature
 * compatibility with the double-precision variants. The unneeded
 * second CUDA vector is discarded, i.e. its memory is not accessed.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::is_same<T, typename cuda_vector_converter<T, tuple<U const, V const> >::result_type>, void>::type
operator<<=(tuple<T&, T&> left, tuple<U&, V&> right)
{
    get<0>(left) = cuda_vector_converter<T, tuple<U const, V const> >::apply(get<0>(right), get<1>(right));
}

/**
 * Pack tuple of double-precision floating-point value and scalar into tuple of two CUDA vectors.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::is_same<tuple<T, T>, typename cuda_vector_converter<tuple<T, T>, tuple<U const, V const> >::result_type>, void>::type
operator<<=(tuple<T&, T&> left, tuple<U&, V&> right)
{
    left = cuda_vector_converter<tuple<T, T>, tuple<U const, V const> >::apply(get<0>(right), get<1>(right));
}

/**
 * Unpack CUDA vector into tuple of single-precision floating-point value and scalar.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::is_same<tuple<U, V>, typename cuda_vector_converter<tuple<U, V>, T const>::result_type>, void>::type
operator<<=(tuple<U&, V&> left, T& right)
{
    left = cuda_vector_converter<tuple<U, V>, T const>::apply(right);
}

/**
 * Unpack CUDA vector into tuple of single-precision floating-point value and scalar.
 *
 * This function accepts a tuple of two CUDA vectors for signature
 * compatibility with the double-precision variants. The unneeded
 * second CUDA vector is discarded, i.e. its memory is not accessed.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::is_same<tuple<U, V>, typename cuda_vector_converter<tuple<U, V>, T const>::result_type>, void>::type
operator<<=(tuple<U&, V&> left, tuple<T&, T&> right)
{
    left = cuda_vector_converter<tuple<U, V>, T const>::apply(get<0>(right));
}

/**
 * Unpack tuple of two CUDA vectors into tuple of double-precision floating-point value and scalar.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::is_same<tuple<U, V>, typename cuda_vector_converter<tuple<U, V>, tuple<T const, T const> >::result_type>, void>::type
operator<<=(tuple<U&, V&> left, tuple<T&, T&> right)
{
    left = cuda_vector_converter<tuple<U, V>, tuple<T const, T const> >::apply(get<0>(right), get<1>(right));
}

/**
 * Pack fixed_vector<float, 3> and int into float4.
 */
template <>
struct cuda_vector_converter<float4, tuple<fixed_vector<float, 3> const, int const> >
{
    typedef float4 result_type;

    static HALMD_GPU_ENABLED result_type apply(fixed_vector<float, 3> const& u, int const& v)
    {
        float4 t;
        t.x = u[0];
        t.y = u[1];
        t.z = u[2];
        t.w = __int_as_float(v);
        return t;
    }
};

/**
 * Unpack float4 into tuple of fixed_vector<float, 3> and int.
 */
template <>
struct cuda_vector_converter<tuple<fixed_vector<float, 3>, int>, float4 const>
{
    typedef tuple<fixed_vector<float, 3>, int> result_type;

    static HALMD_GPU_ENABLED result_type apply(float4 const& t)
    {
        fixed_vector<float, 3> u;
        u[0] = t.x;
        u[1] = t.y;
        u[2] = t.z;
        int v = __float_as_int(t.w);
        return make_tuple(u, v);
    }
};

/**
 * Pack fixed_vector<float, 2> and int into float4.
 */
template <>
struct cuda_vector_converter<float4, tuple<fixed_vector<float, 2> const, int const> >
{
    typedef float4 result_type;

    static HALMD_GPU_ENABLED result_type apply(fixed_vector<float, 2> const& u, int const& v)
    {
        float4 t;
        t.x = u[0];
        t.y = u[1];
        t.z = 0;
        t.w = __int_as_float(v);
        return t;
    }
};

/**
 * Unpack float4 into tuple of fixed_vector<float, 2> and int.
 */
template <>
struct cuda_vector_converter<tuple<fixed_vector<float, 2>, int>, float4 const>
{
    typedef tuple<fixed_vector<float, 2>, int> result_type;

    static HALMD_GPU_ENABLED result_type apply(float4 const& t)
    {
        fixed_vector<float, 2> u;
        u[0] = t.x;
        u[1] = t.y;
        int v = __float_as_int(t.w);
        return make_tuple(u, v);
    }
};

/**
 * Pack fixed_vector<float, 3> and unsigned int into float4.
 */
template <>
struct cuda_vector_converter<float4, tuple<fixed_vector<float, 3> const, unsigned int const> >
{
    typedef float4 result_type;

    static HALMD_GPU_ENABLED result_type apply(fixed_vector<float, 3> const& u, unsigned int const& v)
    {
        float4 t;
        t.x = u[0];
        t.y = u[1];
        t.z = u[2];
        t.w = __int_as_float(v);
        return t;
    }
};

/**
 * Unpack float4 into tuple of fixed_vector<float, 3> and unsigned int.
 */
template <>
struct cuda_vector_converter<tuple<fixed_vector<float, 3>, unsigned int>, float4 const>
{
    typedef tuple<fixed_vector<float, 3>, unsigned int> result_type;

    static HALMD_GPU_ENABLED result_type apply(float4 const& t)
    {
        fixed_vector<float, 3> u;
        u[0] = t.x;
        u[1] = t.y;
        u[2] = t.z;
        unsigned int v = __float_as_int(t.w);
        return make_tuple(u, v);
    }
};

/**
 * Pack fixed_vector<float, 2> and unsigned int into float4.
 */
template <>
struct cuda_vector_converter<float4, tuple<fixed_vector<float, 2> const, unsigned int const> >
{
    typedef float4 result_type;

    static HALMD_GPU_ENABLED result_type apply(fixed_vector<float, 2> const& u, unsigned int const& v)
    {
        float4 t;
        t.x = u[0];
        t.y = u[1];
        t.z = 0;
        t.w = __int_as_float(v);
        return t;
    }
};

/**
 * Unpack float4 into tuple of fixed_vector<float, 2> and unsigned int.
 */
template <>
struct cuda_vector_converter<tuple<fixed_vector<float, 2>, unsigned int>, float4 const>
{
    typedef tuple<fixed_vector<float, 2>, unsigned int> result_type;

    static HALMD_GPU_ENABLED result_type apply(float4 const& t)
    {
        fixed_vector<float, 2> u;
        u[0] = t.x;
        u[1] = t.y;
        unsigned int v = __float_as_int(t.w);
        return make_tuple(u, v);
    }
};

/**
 * Pack fixed_vector<float, 3> and float into float4.
 */
template <>
struct cuda_vector_converter<float4, tuple<fixed_vector<float, 3> const, float const> >
{
    typedef float4 result_type;

    static HALMD_GPU_ENABLED result_type apply(fixed_vector<float, 3> const& u, float const& v)
    {
        float4 t;
        t.x = u[0];
        t.y = u[1];
        t.z = u[2];
        t.w = v;
        return t;
    }
};

/**
 * Unpack float4 into tuple of fixed_vector<float, 3> and float.
 */
template <>
struct cuda_vector_converter<tuple<fixed_vector<float, 3>, float>, float4 const>
{
    typedef tuple<fixed_vector<float, 3>, float> result_type;

    static HALMD_GPU_ENABLED result_type apply(float4 const& t)
    {
        fixed_vector<float, 3> u;
        u[0] = t.x;
        u[1] = t.y;
        u[2] = t.z;
        float v = t.w;
        return make_tuple(u, v);
    }
};

/**
 * Pack fixed_vector<float, 2> and float into float4.
 */
template <>
struct cuda_vector_converter<float4, tuple<fixed_vector<float, 2> const, float const> >
{
    typedef float4 result_type;

    static HALMD_GPU_ENABLED result_type apply(fixed_vector<float, 2> const& u, float const& v)
    {
        float4 t;
        t.x = u[0];
        t.y = u[1];
        t.z = 0;
        t.w = v;
        return t;
    }
};

/**
 * Unpack float4 into tuple of fixed_vector<float, 2> and float.
 */
template <>
struct cuda_vector_converter<tuple<fixed_vector<float, 2>, float>, float4 const>
{
    typedef tuple<fixed_vector<float, 2>, float> result_type;

    static HALMD_GPU_ENABLED result_type apply(float4 const& t)
    {
        fixed_vector<float, 2> u;
        u[0] = t.x;
        u[1] = t.y;
        float v = t.w;
        return make_tuple(u, v);
    }
};

/**
 * Pack fixed_vector<double, 3> and int into tuple of float4 and float4.
 */
template <>
struct cuda_vector_converter<tuple<float4, float4>, tuple<fixed_vector<double, 3> const, int const> >
{
    typedef tuple<float4, float4> result_type;

    static HALMD_GPU_ENABLED result_type apply(fixed_vector<double, 3> const& u, int const& v)
    {
        float4 hi, lo;
        tie(hi.x, lo.x) = split(dsfloat(u[0]));
        tie(hi.y, lo.y) = split(dsfloat(u[1]));
        tie(hi.z, lo.z) = split(dsfloat(u[2]));
        tie(hi.w, lo.w) = make_tuple(__int_as_float(v), 0);
        return make_tuple(hi, lo);
    }
};

/**
 * Unpack float4 and float4 into tuple of fixed_vector<double, 3> and int.
 */
template <>
struct cuda_vector_converter<tuple<fixed_vector<double, 3>, int>, tuple<float4 const, float4 const> >
{
    typedef tuple<fixed_vector<double, 3>, int> result_type;

    static HALMD_GPU_ENABLED result_type apply(float4 const& hi, float4 const& lo)
    {
        fixed_vector<double, 3> u;
        u[0] = dsfloat(hi.x, lo.x);
        u[1] = dsfloat(hi.y, lo.y);
        u[2] = dsfloat(hi.z, lo.z);
        int v = __float_as_int(hi.w);
        return make_tuple(u, v);
    }
};

/**
 * Pack fixed_vector<double, 2> and int into tuple of float4 and float4.
 */
template <>
struct cuda_vector_converter<tuple<float4, float4>, tuple<fixed_vector<double, 2> const, int const> >
{
    typedef tuple<float4, float4> result_type;

    static HALMD_GPU_ENABLED result_type apply(fixed_vector<double, 2> const& u, int const& v)
    {
        float4 hi, lo;
        tie(hi.x, lo.x) = split(dsfloat(u[0]));
        tie(hi.y, lo.y) = split(dsfloat(u[1]));
        tie(hi.z, lo.z) = make_tuple(0, 0);
        tie(hi.w, lo.w) = make_tuple(__int_as_float(v), 0);
        return make_tuple(hi, lo);
    }
};

/**
 * Unpack float4 and float4 into tuple of fixed_vector<double, 2> and int.
 */
template <>
struct cuda_vector_converter<tuple<fixed_vector<double, 2>, int>, tuple<float4 const, float4 const> >
{
    typedef tuple<fixed_vector<double, 2>, int> result_type;

    static HALMD_GPU_ENABLED result_type apply(float4 const& hi, float4 const& lo)
    {
        fixed_vector<double, 2> u;
        u[0] = dsfloat(hi.x, lo.x);
        u[1] = dsfloat(hi.y, lo.y);
        int v = __float_as_int(hi.w);
        return make_tuple(u, v);
    }
};

/**
 * Pack fixed_vector<double, 3> and unsigned int into tuple of float4 and float4.
 */
template <>
struct cuda_vector_converter<tuple<float4, float4>, tuple<fixed_vector<double, 3> const, unsigned int const> >
{
    typedef tuple<float4, float4> result_type;

    static HALMD_GPU_ENABLED result_type apply(fixed_vector<double, 3> const& u, unsigned int const& v)
    {
        float4 hi, lo;
        tie(hi.x, lo.x) = split(dsfloat(u[0]));
        tie(hi.y, lo.y) = split(dsfloat(u[1]));
        tie(hi.z, lo.z) = split(dsfloat(u[2]));
        tie(hi.w, lo.w) = make_tuple(__int_as_float(v), 0);
        return make_tuple(hi, lo);
    }
};

/**
 * Unpack float4 and float4 into tuple of fixed_vector<double, 3> and unsigned int.
 */
template <>
struct cuda_vector_converter<tuple<fixed_vector<double, 3>, unsigned int>, tuple<float4 const, float4 const> >
{
    typedef tuple<fixed_vector<double, 3>, unsigned int> result_type;

    static HALMD_GPU_ENABLED result_type apply(float4 const& hi, float4 const& lo)
    {
        fixed_vector<double, 3> u;
        u[0] = dsfloat(hi.x, lo.x);
        u[1] = dsfloat(hi.y, lo.y);
        u[2] = dsfloat(hi.z, lo.z);
        unsigned int v = __float_as_int(hi.w);
        return make_tuple(u, v);
    }
};

/**
 * Pack fixed_vector<double, 2> and unsigned int into tuple of float4 and float4.
 */
template <>
struct cuda_vector_converter<tuple<float4, float4>, tuple<fixed_vector<double, 2> const, unsigned int const> >
{
    typedef tuple<float4, float4> result_type;

    static HALMD_GPU_ENABLED result_type apply(fixed_vector<double, 2> const& u, unsigned int const& v)
    {
        float4 hi, lo;
        tie(hi.x, lo.x) = split(dsfloat(u[0]));
        tie(hi.y, lo.y) = split(dsfloat(u[1]));
        tie(hi.z, lo.z) = make_tuple(0, 0);
        tie(hi.w, lo.w) = make_tuple(__int_as_float(v), 0);
        return make_tuple(hi, lo);
    }
};

/**
 * Unpack float4 and float4 into tuple of fixed_vector<double, 2> and unsigned int.
 */
template <>
struct cuda_vector_converter<tuple<fixed_vector<double, 2>, unsigned int>, tuple<float4 const, float4 const> >
{
    typedef tuple<fixed_vector<double, 2>, unsigned int> result_type;

    static HALMD_GPU_ENABLED result_type apply(float4 const& hi, float4 const& lo)
    {
        fixed_vector<double, 2> u;
        u[0] = dsfloat(hi.x, lo.x);
        u[1] = dsfloat(hi.y, lo.y);
        unsigned int v = __float_as_int(hi.w);
        return make_tuple(u, v);
    }
};

/**
 * Pack fixed_vector<double, 3> and float into tuple of float4 and float4.
 */
template <>
struct cuda_vector_converter<tuple<float4, float4>, tuple<fixed_vector<double, 3> const, float const> >
{
    typedef tuple<float4, float4> result_type;

    static HALMD_GPU_ENABLED result_type apply(fixed_vector<double, 3> const& u, float const& v)
    {
        float4 hi, lo;
        tie(hi.x, lo.x) = split(dsfloat(u[0]));
        tie(hi.y, lo.y) = split(dsfloat(u[1]));
        tie(hi.z, lo.z) = split(dsfloat(u[2]));
        tie(hi.w, lo.w) = make_tuple(v, 0);
        return make_tuple(hi, lo);
    }
};

/**
 * Unpack float4 and float4 into tuple of fixed_vector<double, 3> and float.
 */
template <>
struct cuda_vector_converter<tuple<fixed_vector<double, 3>, float>, tuple<float4 const, float4 const> >
{
    typedef tuple<fixed_vector<double, 3>, float> result_type;

    static HALMD_GPU_ENABLED result_type apply(float4 const& hi, float4 const& lo)
    {
        fixed_vector<double, 3> u;
        u[0] = dsfloat(hi.x, lo.x);
        u[1] = dsfloat(hi.y, lo.y);
        u[2] = dsfloat(hi.z, lo.z);
        float v = hi.w;
        return make_tuple(u, v);
    }
};

/**
 * Pack fixed_vector<double, 2> and float into tuple of float4 and float4.
 */
template <>
struct cuda_vector_converter<tuple<float4, float4>, tuple<fixed_vector<double, 2> const, float const> >
{
    typedef tuple<float4, float4> result_type;

    static HALMD_GPU_ENABLED result_type apply(fixed_vector<double, 2> const& u, float const& v)
    {
        float4 hi, lo;
        tie(hi.x, lo.x) = split(dsfloat(u[0]));
        tie(hi.y, lo.y) = split(dsfloat(u[1]));
        tie(hi.z, lo.z) = make_tuple(0, 0);
        tie(hi.w, lo.w) = make_tuple(v, 0);
        return make_tuple(hi, lo);
    }
};

/**
 * Unpack float4 and float4 into tuple of fixed_vector<double, 2> and float.
 */
template <>
struct cuda_vector_converter<tuple<fixed_vector<double, 2>, float>, tuple<float4 const, float4 const> >
{
    typedef tuple<fixed_vector<double, 2>, float> result_type;

    static HALMD_GPU_ENABLED result_type apply(float4 const& hi, float4 const& lo)
    {
        fixed_vector<double, 2> u;
        u[0] = dsfloat(hi.x, lo.x);
        u[1] = dsfloat(hi.y, lo.y);
        float v = hi.w;
        return make_tuple(u, v);
    }
};

/**
 * Pack fixed_vector<double, 3> and double into tuple of float4 and float4.
 */
template <>
struct cuda_vector_converter<tuple<float4, float4>, tuple<fixed_vector<double, 3> const, double const> >
{
    typedef tuple<float4, float4> result_type;

    static HALMD_GPU_ENABLED result_type apply(fixed_vector<double, 3> const& u, double const& v)
    {
        float4 hi, lo;
        tie(hi.x, lo.x) = split(dsfloat(u[0]));
        tie(hi.y, lo.y) = split(dsfloat(u[1]));
        tie(hi.z, lo.z) = split(dsfloat(u[2]));
        tie(hi.w, lo.w) = split(dsfloat(v));
        return make_tuple(hi, lo);
    }
};

/**
 * Unpack float4 and float4 into tuple of fixed_vector<double, 3> and double.
 */
template <>
struct cuda_vector_converter<tuple<fixed_vector<double, 3>, double>, tuple<float4 const, float4 const> >
{
    typedef tuple<fixed_vector<double, 3>, double> result_type;

    static HALMD_GPU_ENABLED result_type apply(float4 const& hi, float4 const& lo)
    {
        fixed_vector<double, 3> u;
        u[0] = dsfloat(hi.x, lo.x);
        u[1] = dsfloat(hi.y, lo.y);
        u[2] = dsfloat(hi.z, lo.z);
        double v = dsfloat(hi.w, lo.w);
        return make_tuple(u, v);
    }
};

/**
 * Pack fixed_vector<double, 2> and double into tuple of float4 and float4.
 */
template <>
struct cuda_vector_converter<tuple<float4, float4>, tuple<fixed_vector<double, 2> const, double const> >
{
    typedef tuple<float4, float4> result_type;

    static HALMD_GPU_ENABLED result_type apply(fixed_vector<double, 2> const& u, double const& v)
    {
        float4 hi, lo;
        tie(hi.x, lo.x) = split(dsfloat(u[0]));
        tie(hi.y, lo.y) = split(dsfloat(u[1]));
        tie(hi.z, lo.z) = make_tuple(0, 0);
        tie(hi.w, lo.w) = split(dsfloat(v));
        return make_tuple(hi, lo);
    }
};

/**
 * Unpack float4 and float4 into tuple of fixed_vector<double, 2> and double.
 */
template <>
struct cuda_vector_converter<tuple<fixed_vector<double, 2>, double>, tuple<float4 const, float4 const> >
{
    typedef tuple<fixed_vector<double, 2>, double> result_type;

    static HALMD_GPU_ENABLED result_type apply(float4 const& hi, float4 const& lo)
    {
        fixed_vector<double, 2> u;
        u[0] = dsfloat(hi.x, lo.x);
        u[1] = dsfloat(hi.y, lo.y);
        double v = dsfloat(hi.w, lo.w);
        return make_tuple(u, v);
    }
};

/**
 * Pack fixed_vector<dsfloat, 3> and int into tuple of float4 and float4.
 */
template <>
struct cuda_vector_converter<tuple<float4, float4>, tuple<fixed_vector<dsfloat, 3> const, int const> >
{
    typedef tuple<float4, float4> result_type;

    static HALMD_GPU_ENABLED result_type apply(fixed_vector<dsfloat, 3> const& u, int const& v)
    {
        float4 hi, lo;
        tie(hi.x, lo.x) = split(u[0]);
        tie(hi.y, lo.y) = split(u[1]);
        tie(hi.z, lo.z) = split(u[2]);
        tie(hi.w, lo.w) = make_tuple(__int_as_float(v), 0);
        return make_tuple(hi, lo);
    }
};

/**
 * Unpack float4 and float4 into tuple of fixed_vector<dsfloat, 3> and int.
 */
template <>
struct cuda_vector_converter<tuple<fixed_vector<dsfloat, 3>, int>, tuple<float4 const, float4 const> >
{
    typedef tuple<fixed_vector<dsfloat, 3>, int> result_type;

    static HALMD_GPU_ENABLED result_type apply(float4 const& hi, float4 const& lo)
    {
        fixed_vector<dsfloat, 3> u;
        u[0] = dsfloat(hi.x, lo.x);
        u[1] = dsfloat(hi.y, lo.y);
        u[2] = dsfloat(hi.z, lo.z);
        int v = __float_as_int(hi.w);
        return make_tuple(u, v);
    }
};

/**
 * Pack fixed_vector<dsfloat, 2> and int into tuple of float4 and float4.
 */
template <>
struct cuda_vector_converter<tuple<float4, float4>, tuple<fixed_vector<dsfloat, 2> const, int const> >
{
    typedef tuple<float4, float4> result_type;

    static HALMD_GPU_ENABLED result_type apply(fixed_vector<dsfloat, 2> const& u, int const& v)
    {
        float4 hi, lo;
        tie(hi.x, lo.x) = split(u[0]);
        tie(hi.y, lo.y) = split(u[1]);
        tie(hi.z, lo.z) = make_tuple(0, 0);
        tie(hi.w, lo.w) = make_tuple(__int_as_float(v), 0);
        return make_tuple(hi, lo);
    }
};

/**
 * Unpack float4 and float4 into tuple of fixed_vector<dsfloat, 2> and int.
 */
template <>
struct cuda_vector_converter<tuple<fixed_vector<dsfloat, 2>, int>, tuple<float4 const, float4 const> >
{
    typedef tuple<fixed_vector<dsfloat, 2>, int> result_type;

    static HALMD_GPU_ENABLED result_type apply(float4 const& hi, float4 const& lo)
    {
        fixed_vector<dsfloat, 2> u;
        u[0] = dsfloat(hi.x, lo.x);
        u[1] = dsfloat(hi.y, lo.y);
        int v = __float_as_int(hi.w);
        return make_tuple(u, v);
    }
};

/**
 * Pack fixed_vector<dsfloat, 3> and unsigned int into tuple of float4 and float4.
 */
template <>
struct cuda_vector_converter<tuple<float4, float4>, tuple<fixed_vector<dsfloat, 3> const, unsigned int const> >
{
    typedef tuple<float4, float4> result_type;

    static HALMD_GPU_ENABLED result_type apply(fixed_vector<dsfloat, 3> const& u, unsigned int const& v)
    {
        float4 hi, lo;
        tie(hi.x, lo.x) = split(u[0]);
        tie(hi.y, lo.y) = split(u[1]);
        tie(hi.z, lo.z) = split(u[2]);
        tie(hi.w, lo.w) = make_tuple(__int_as_float(v), 0);
        return make_tuple(hi, lo);
    }
};

/**
 * Unpack float4 and float4 into tuple of fixed_vector<dsfloat, 3> and unsigned int.
 */
template <>
struct cuda_vector_converter<tuple<fixed_vector<dsfloat, 3>, unsigned int>, tuple<float4 const, float4 const> >
{
    typedef tuple<fixed_vector<dsfloat, 3>, unsigned int> result_type;

    static HALMD_GPU_ENABLED result_type apply(float4 const& hi, float4 const& lo)
    {
        fixed_vector<dsfloat, 3> u;
        u[0] = dsfloat(hi.x, lo.x);
        u[1] = dsfloat(hi.y, lo.y);
        u[2] = dsfloat(hi.z, lo.z);
        unsigned int v = __float_as_int(hi.w);
        return make_tuple(u, v);
    }
};

/**
 * Pack fixed_vector<dsfloat, 2> and unsigned int into tuple of float4 and float4.
 */
template <>
struct cuda_vector_converter<tuple<float4, float4>, tuple<fixed_vector<dsfloat, 2> const, unsigned int const> >
{
    typedef tuple<float4, float4> result_type;

    static HALMD_GPU_ENABLED result_type apply(fixed_vector<dsfloat, 2> const& u, unsigned int const& v)
    {
        float4 hi, lo;
        tie(hi.x, lo.x) = split(u[0]);
        tie(hi.y, lo.y) = split(u[1]);
        tie(hi.z, lo.z) = make_tuple(0, 0);
        tie(hi.w, lo.w) = make_tuple(__int_as_float(v), 0);
        return make_tuple(hi, lo);
    }
};

/**
 * Unpack float4 and float4 into tuple of fixed_vector<dsfloat, 2> and unsigned int.
 */
template <>
struct cuda_vector_converter<tuple<fixed_vector<dsfloat, 2>, unsigned int>, tuple<float4 const, float4 const> >
{
    typedef tuple<fixed_vector<dsfloat, 2>, unsigned int> result_type;

    static HALMD_GPU_ENABLED result_type apply(float4 const& hi, float4 const& lo)
    {
        fixed_vector<dsfloat, 2> u;
        u[0] = dsfloat(hi.x, lo.x);
        u[1] = dsfloat(hi.y, lo.y);
        unsigned int v = __float_as_int(hi.w);
        return make_tuple(u, v);
    }
};

/**
 * Pack fixed_vector<dsfloat, 3> and float into tuple of float4 and float4.
 */
template <>
struct cuda_vector_converter<tuple<float4, float4>, tuple<fixed_vector<dsfloat, 3> const, float const> >
{
    typedef tuple<float4, float4> result_type;

    static HALMD_GPU_ENABLED result_type apply(fixed_vector<dsfloat, 3> const& u, float const& v)
    {
        float4 hi, lo;
        tie(hi.x, lo.x) = split(u[0]);
        tie(hi.y, lo.y) = split(u[1]);
        tie(hi.z, lo.z) = split(u[2]);
        tie(hi.w, lo.w) = make_tuple(v, 0);
        return make_tuple(hi, lo);
    }
};

/**
 * Unpack float4 and float4 into tuple of fixed_vector<dsfloat, 3> and float.
 */
template <>
struct cuda_vector_converter<tuple<fixed_vector<dsfloat, 3>, float>, tuple<float4 const, float4 const> >
{
    typedef tuple<fixed_vector<dsfloat, 3>, float> result_type;

    static HALMD_GPU_ENABLED result_type apply(float4 const& hi, float4 const& lo)
    {
        fixed_vector<dsfloat, 3> u;
        u[0] = dsfloat(hi.x, lo.x);
        u[1] = dsfloat(hi.y, lo.y);
        u[2] = dsfloat(hi.z, lo.z);
        float v = hi.w;
        return make_tuple(u, v);
    }
};

/**
 * Pack fixed_vector<dsfloat, 2> and float into tuple of float4 and float4.
 */
template <>
struct cuda_vector_converter<tuple<float4, float4>, tuple<fixed_vector<dsfloat, 2> const, float const> >
{
    typedef tuple<float4, float4> result_type;

    static HALMD_GPU_ENABLED result_type apply(fixed_vector<dsfloat, 2> const& u, float const& v)
    {
        float4 hi, lo;
        tie(hi.x, lo.x) = split(u[0]);
        tie(hi.y, lo.y) = split(u[1]);
        tie(hi.z, lo.z) = make_tuple(0, 0);
        tie(hi.w, lo.w) = make_tuple(v, 0);
        return make_tuple(hi, lo);
    }
};

/**
 * Unpack float4 and float4 into tuple of fixed_vector<dsfloat, 2> and float.
 */
template <>
struct cuda_vector_converter<tuple<fixed_vector<dsfloat, 2>, float>, tuple<float4 const, float4 const> >
{
    typedef tuple<fixed_vector<dsfloat, 2>, float> result_type;

    static HALMD_GPU_ENABLED result_type apply(float4 const& hi, float4 const& lo)
    {
        fixed_vector<dsfloat, 2> u;
        u[0] = dsfloat(hi.x, lo.x);
        u[1] = dsfloat(hi.y, lo.y);
        float v = hi.w;
        return make_tuple(u, v);
    }
};

/**
 * Pack fixed_vector<dsfloat, 3> and dsfloat into tuple of float4 and float4.
 */
template <>
struct cuda_vector_converter<tuple<float4, float4>, tuple<fixed_vector<dsfloat, 3> const, dsfloat const> >
{
    typedef tuple<float4, float4> result_type;

    static HALMD_GPU_ENABLED result_type apply(fixed_vector<dsfloat, 3> const& u, dsfloat const& v)
    {
        float4 hi, lo;
        tie(hi.x, lo.x) = split(u[0]);
        tie(hi.y, lo.y) = split(u[1]);
        tie(hi.z, lo.z) = split(u[2]);
        tie(hi.w, lo.w) = split(v);
        return make_tuple(hi, lo);
    }
};

/**
 * Unpack float4 and float4 into tuple of fixed_vector<dsfloat, 3> and dsfloat.
 */
template <>
struct cuda_vector_converter<tuple<fixed_vector<dsfloat, 3>, dsfloat>, tuple<float4 const, float4 const> >
{
    typedef tuple<fixed_vector<dsfloat, 3>, dsfloat> result_type;

    static HALMD_GPU_ENABLED result_type apply(float4 const& hi, float4 const& lo)
    {
        fixed_vector<dsfloat, 3> u;
        u[0] = dsfloat(hi.x, lo.x);
        u[1] = dsfloat(hi.y, lo.y);
        u[2] = dsfloat(hi.z, lo.z);
        dsfloat v = dsfloat(hi.w, lo.w);
        return make_tuple(u, v);
    }
};

/**
 * Pack fixed_vector<dsfloat, 2> and dsfloat into tuple of float4 and float4.
 */
template <>
struct cuda_vector_converter<tuple<float4, float4>, tuple<fixed_vector<dsfloat, 2> const, dsfloat const> >
{
    typedef tuple<float4, float4> result_type;

    static HALMD_GPU_ENABLED result_type apply(fixed_vector<dsfloat, 2> const& u, dsfloat const& v)
    {
        float4 hi, lo;
        tie(hi.x, lo.x) = split(u[0]);
        tie(hi.y, lo.y) = split(u[1]);
        tie(hi.z, lo.z) = make_tuple(0, 0);
        tie(hi.w, lo.w) = split(v);
        return make_tuple(hi, lo);
    }
};

/**
 * Unpack float4 and float4 into tuple of fixed_vector<dsfloat, 2> and dsfloat.
 */
template <>
struct cuda_vector_converter<tuple<fixed_vector<dsfloat, 2>, dsfloat>, tuple<float4 const, float4 const> >
{
    typedef tuple<fixed_vector<dsfloat, 2>, dsfloat> result_type;

    static HALMD_GPU_ENABLED result_type apply(float4 const& hi, float4 const& lo)
    {
        fixed_vector<dsfloat, 2> u;
        u[0] = dsfloat(hi.x, lo.x);
        u[1] = dsfloat(hi.y, lo.y);
        dsfloat v = dsfloat(hi.w, lo.w);
        return make_tuple(u, v);
    }
};

} // namespace detail
} // namespace numeric
} // namespace blas
} // namespace halmd

#endif /* ! HALMD_NUMERIC_BLAS_FIXED_VECTOR_CUDA_VECTOR_CONVERTER_HPP */
