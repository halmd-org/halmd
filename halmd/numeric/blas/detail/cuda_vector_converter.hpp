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

#ifndef HALMD_NUMERIC_BLAS_DETAIL_CUDA_VECTOR_CONVERTER_HPP
#define HALMD_NUMERIC_BLAS_DETAIL_CUDA_VECTOR_CONVERTER_HPP

#include <halmd/config.hpp>

#include <boost/mpl/and.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>

#include <halmd/numeric/blas/detail/vector.hpp>
#include <halmd/numeric/mp/dsfloat.hpp>
#include <halmd/utility/tuple.hpp>

namespace halmd {
namespace numeric {
namespace blas {
namespace detail {

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

/**
 * Pack fixed_vector<float, 3> and int into float4.
 *
 * Takes right-hand tuple of references, as returned by tie.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U const, fixed_vector<float, 3> const>
      , boost::is_same<V const, int const>
    >, void>::type
operator<<=(T& left, tuple<U&, V&> right)
{
    fixed_vector<float, 3> u;
    int v;
    tie(u, v) = right;
    float4 t;
    t.x = u[0];
    t.y = u[1];
    t.z = u[2];
    t.w = __int_as_float(v);
    left = t;
}

/**
 * Pack fixed_vector<float, 3> and int into float4.
 *
 * Takes right-hand tuple of values, as returned by make_tuple.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<float, 3> >
      , boost::is_same<V, int>
    >, void>::type
operator<<=(T& left, tuple<U, V> right)
{
    fixed_vector<float, 3> u;
    int v;
    tie(u, v) = right;
    float4 t;
    t.x = u[0];
    t.y = u[1];
    t.z = u[2];
    t.w = __int_as_float(v);
    left = t;
}

/**
 * Unpack float4 into tuple of fixed_vector<float, 3> and int.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<float, 3> >
      , boost::is_same<V, int>
    >, void>::type
operator<<=(tuple<U&, V&> left, T const& right)
{
    fixed_vector<float, 3> u;
    u[0] = right.x;
    u[1] = right.y;
    u[2] = right.z;
    int v = __float_as_int(right.w);
    left = make_tuple(u, v);
}

/**
 * Pack fixed_vector<float, 2> and int into float4.
 *
 * Takes right-hand tuple of references, as returned by tie.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U const, fixed_vector<float, 2> const>
      , boost::is_same<V const, int const>
    >, void>::type
operator<<=(T& left, tuple<U&, V&> right)
{
    fixed_vector<float, 2> u;
    int v;
    tie(u, v) = right;
    float4 t;
    t.x = u[0];
    t.y = u[1];
    t.z = 0;
    t.w = __int_as_float(v);
    left = t;
}

/**
 * Pack fixed_vector<float, 2> and int into float4.
 *
 * Takes right-hand tuple of values, as returned by make_tuple.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<float, 2> >
      , boost::is_same<V, int>
    >, void>::type
operator<<=(T& left, tuple<U, V> right)
{
    fixed_vector<float, 2> u;
    int v;
    tie(u, v) = right;
    float4 t;
    t.x = u[0];
    t.y = u[1];
    t.z = 0;
    t.w = __int_as_float(v);
    left = t;
}

/**
 * Unpack float4 into tuple of fixed_vector<float, 2> and int.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<float, 2> >
      , boost::is_same<V, int>
    >, void>::type
operator<<=(tuple<U&, V&> left, T const& right)
{
    fixed_vector<float, 2> u;
    u[0] = right.x;
    u[1] = right.y;
    int v = __float_as_int(right.w);
    left = make_tuple(u, v);
}

/**
 * Pack fixed_vector<float, 3> and unsigned int into float4.
 *
 * Takes right-hand tuple of references, as returned by tie.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U const, fixed_vector<float, 3> const>
      , boost::is_same<V const, unsigned int const>
    >, void>::type
operator<<=(T& left, tuple<U&, V&> right)
{
    fixed_vector<float, 3> u;
    unsigned int v;
    tie(u, v) = right;
    float4 t;
    t.x = u[0];
    t.y = u[1];
    t.z = u[2];
    t.w = __int_as_float(v);
    left = t;
}

/**
 * Pack fixed_vector<float, 3> and unsigned int into float4.
 *
 * Takes right-hand tuple of values, as returned by make_tuple.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<float, 3> >
      , boost::is_same<V, unsigned int>
    >, void>::type
operator<<=(T& left, tuple<U, V> right)
{
    fixed_vector<float, 3> u;
    unsigned int v;
    tie(u, v) = right;
    float4 t;
    t.x = u[0];
    t.y = u[1];
    t.z = u[2];
    t.w = __int_as_float(v);
    left = t;
}

/**
 * Unpack float4 into tuple of fixed_vector<float, 3> and unsigned int.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<float, 3> >
      , boost::is_same<V, unsigned int>
    >, void>::type
operator<<=(tuple<U&, V&> left, T const& right)
{
    fixed_vector<float, 3> u;
    u[0] = right.x;
    u[1] = right.y;
    u[2] = right.z;
    unsigned int v = __float_as_int(right.w);
    left = make_tuple(u, v);
}

/**
 * Pack fixed_vector<float, 2> and unsigned int into float4.
 *
 * Takes right-hand tuple of references, as returned by tie.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U const, fixed_vector<float, 2> const>
      , boost::is_same<V const, unsigned int const>
    >, void>::type
operator<<=(T& left, tuple<U&, V&> right)
{
    fixed_vector<float, 2> u;
    unsigned int v;
    tie(u, v) = right;
    float4 t;
    t.x = u[0];
    t.y = u[1];
    t.z = 0;
    t.w = __int_as_float(v);
    left = t;
}

/**
 * Pack fixed_vector<float, 2> and unsigned int into float4.
 *
 * Takes right-hand tuple of values, as returned by make_tuple.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<float, 2> >
      , boost::is_same<V, unsigned int>
    >, void>::type
operator<<=(T& left, tuple<U, V> right)
{
    fixed_vector<float, 2> u;
    unsigned int v;
    tie(u, v) = right;
    float4 t;
    t.x = u[0];
    t.y = u[1];
    t.z = 0;
    t.w = __int_as_float(v);
    left = t;
}

/**
 * Unpack float4 into tuple of fixed_vector<float, 2> and unsigned int.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<float, 2> >
      , boost::is_same<V, unsigned int>
    >, void>::type
operator<<=(tuple<U&, V&> left, T const& right)
{
    fixed_vector<float, 2> u;
    u[0] = right.x;
    u[1] = right.y;
    unsigned int v = __float_as_int(right.w);
    left = make_tuple(u, v);
}

/**
 * Pack fixed_vector<float, 3> and float into float4.
 *
 * Takes right-hand tuple of references, as returned by tie.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U const, fixed_vector<float, 3> const>
      , boost::is_same<V const, float const>
    >, void>::type
operator<<=(T& left, tuple<U&, V&> right)
{
    fixed_vector<float, 3> u;
    float v;
    tie(u, v) = right;
    float4 t;
    t.x = u[0];
    t.y = u[1];
    t.z = u[2];
    t.w = v;
    left = t;
}

/**
 * Pack fixed_vector<float, 3> and float into float4.
 *
 * Takes right-hand tuple of values, as returned by make_tuple.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<float, 3> >
      , boost::is_same<V, float>
    >, void>::type
operator<<=(T& left, tuple<U, V> right)
{
    fixed_vector<float, 3> u;
    float v;
    tie(u, v) = right;
    float4 t;
    t.x = u[0];
    t.y = u[1];
    t.z = u[2];
    t.w = v;
    left = t;
}

/**
 * Unpack float4 into tuple of fixed_vector<float, 3> and float.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<float, 3> >
      , boost::is_same<V, float>
    >, void>::type
operator<<=(tuple<U&, V&> left, T const& right)
{
    fixed_vector<float, 3> u;
    u[0] = right.x;
    u[1] = right.y;
    u[2] = right.z;
    float v = right.w;
    left = make_tuple(u, v);
}

/**
 * Pack fixed_vector<float, 2> and float into float4.
 *
 * Takes right-hand tuple of references, as returned by tie.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U const, fixed_vector<float, 2> const>
      , boost::is_same<V const, float const>
    >, void>::type
operator<<=(T& left, tuple<U&, V&> right)
{
    fixed_vector<float, 2> u;
    float v;
    tie(u, v) = right;
    float4 t;
    t.x = u[0];
    t.y = u[1];
    t.z = 0;
    t.w = v;
    left = t;
}

/**
 * Pack fixed_vector<float, 2> and float into float4.
 *
 * Takes right-hand tuple of values, as returned by make_tuple.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<float, 2> >
      , boost::is_same<V, float>
    >, void>::type
operator<<=(T& left, tuple<U, V> right)
{
    fixed_vector<float, 2> u;
    float v;
    tie(u, v) = right;
    float4 t;
    t.x = u[0];
    t.y = u[1];
    t.z = 0;
    t.w = v;
    left = t;
}

/**
 * Unpack float4 into tuple of fixed_vector<float, 2> and float.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<float, 2> >
      , boost::is_same<V, float>
    >, void>::type
operator<<=(tuple<U&, V&> left, T const& right)
{
    fixed_vector<float, 2> u;
    u[0] = right.x;
    u[1] = right.y;
    float v = right.w;
    left = make_tuple(u, v);
}

/**
 * Pack fixed_vector<double, 3> and int into tuple of float4 and float4.
 *
 * Takes right-hand tuple of references, as returned by tie.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U const, fixed_vector<double, 3> const>
      , boost::is_same<V const, int const>
    >, void>::type
operator<<=(tuple<T&, T&> left, tuple<U&, V&> right)
{
    fixed_vector<double, 3> u;
    int v;
    tie(u, v) = right;
    float4 hi, lo;
    tie(hi.x, lo.x) = split(dsfloat(u[0]));
    tie(hi.y, lo.y) = split(dsfloat(u[1]));
    tie(hi.z, lo.z) = split(dsfloat(u[2]));
    tie(hi.w, lo.w) = make_tuple(__int_as_float(v), 0);
    left = make_tuple(hi, lo);
}

/**
 * Pack fixed_vector<double, 3> and int into tuple of float4 and float4.
 *
 * Takes right-hand tuple of values, as returned by make_tuple.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<double, 3> >
      , boost::is_same<V, int>
    >, void>::type
operator<<=(tuple<T&, T&> left, tuple<U, V> right)
{
    fixed_vector<double, 3> u;
    int v;
    tie(u, v) = right;
    float4 hi, lo;
    tie(hi.x, lo.x) = split(dsfloat(u[0]));
    tie(hi.y, lo.y) = split(dsfloat(u[1]));
    tie(hi.z, lo.z) = split(dsfloat(u[2]));
    tie(hi.w, lo.w) = make_tuple(__int_as_float(v), 0);
    left = make_tuple(hi, lo);
}

/**
 * Unpack float4 and float4 into tuple of fixed_vector<double, 3> and int.
 *
 * Takes right-hand tuple of references, as returned by tie.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T const, float4 const>
      , boost::is_same<U, fixed_vector<double, 3> >
      , boost::is_same<V, int>
    >, void>::type
operator<<=(tuple<U&, V&> left, tuple<T&, T&> right)
{
    float4 hi, lo;
    tie(hi, lo) = right;
    fixed_vector<double, 3> u;
    u[0] = dsfloat(hi.x, lo.x);
    u[1] = dsfloat(hi.y, lo.y);
    u[2] = dsfloat(hi.z, lo.z);
    int v = __float_as_int(hi.w);
    left = make_tuple(u, v);
}

/**
 * Unpack float4 and float4 into tuple of fixed_vector<double, 3> and int.
 *
 * Takes right-hand tuple of values, as returned by make_tuple.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<double, 3> >
      , boost::is_same<V, int>
    >, void>::type
operator<<=(tuple<U&, V&> left, tuple<T, T> right)
{
    float4 hi, lo;
    tie(hi, lo) = right;
    fixed_vector<double, 3> u;
    u[0] = dsfloat(hi.x, lo.x);
    u[1] = dsfloat(hi.y, lo.y);
    u[2] = dsfloat(hi.z, lo.z);
    int v = __float_as_int(hi.w);
    left = make_tuple(u, v);
}

/**
 * Pack fixed_vector<double, 2> and int into tuple of float4 and float4.
 *
 * Takes right-hand tuple of references, as returned by tie.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U const, fixed_vector<double, 2> const>
      , boost::is_same<V const, int const>
    >, void>::type
operator<<=(tuple<T&, T&> left, tuple<U&, V&> right)
{
    fixed_vector<double, 2> u;
    int v;
    tie(u, v) = right;
    float4 hi, lo;
    tie(hi.x, lo.x) = split(dsfloat(u[0]));
    tie(hi.y, lo.y) = split(dsfloat(u[1]));
    tie(hi.z, lo.z) = make_tuple(0, 0);
    tie(hi.w, lo.w) = make_tuple(__int_as_float(v), 0);
    left = make_tuple(hi, lo);
}

/**
 * Pack fixed_vector<double, 2> and int into tuple of float4 and float4.
 *
 * Takes right-hand tuple of values, as returned by make_tuple.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<double, 2> >
      , boost::is_same<V, int>
    >, void>::type
operator<<=(tuple<T&, T&> left, tuple<U, V> right)
{
    fixed_vector<double, 2> u;
    int v;
    tie(u, v) = right;
    float4 hi, lo;
    tie(hi.x, lo.x) = split(dsfloat(u[0]));
    tie(hi.y, lo.y) = split(dsfloat(u[1]));
    tie(hi.z, lo.z) = make_tuple(0, 0);
    tie(hi.w, lo.w) = make_tuple(__int_as_float(v), 0);
    left = make_tuple(hi, lo);
}

/**
 * Unpack float4 and float4 into tuple of fixed_vector<double, 2> and int.
 *
 * Takes right-hand tuple of references, as returned by tie.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T const, float4 const>
      , boost::is_same<U, fixed_vector<double, 2> >
      , boost::is_same<V, int>
    >, void>::type
operator<<=(tuple<U&, V&> left, tuple<T&, T&> right)
{
    float4 hi, lo;
    tie(hi, lo) = right;
    fixed_vector<double, 2> u;
    u[0] = dsfloat(hi.x, lo.x);
    u[1] = dsfloat(hi.y, lo.y);
    int v = __float_as_int(hi.w);
    left = make_tuple(u, v);
}

/**
 * Unpack float4 and float4 into tuple of fixed_vector<double, 2> and int.
 *
 * Takes right-hand tuple of values, as returned by make_tuple.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<double, 2> >
      , boost::is_same<V, int>
    >, void>::type
operator<<=(tuple<U&, V&> left, tuple<T, T> right)
{
    float4 hi, lo;
    tie(hi, lo) = right;
    fixed_vector<double, 2> u;
    u[0] = dsfloat(hi.x, lo.x);
    u[1] = dsfloat(hi.y, lo.y);
    int v = __float_as_int(hi.w);
    left = make_tuple(u, v);
}

/**
 * Pack fixed_vector<double, 3> and unsigned int into tuple of float4 and float4.
 *
 * Takes right-hand tuple of references, as returned by tie.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U const, fixed_vector<double, 3> const>
      , boost::is_same<V const, unsigned int const>
    >, void>::type
operator<<=(tuple<T&, T&> left, tuple<U&, V&> right)
{
    fixed_vector<double, 3> u;
    unsigned int v;
    tie(u, v) = right;
    float4 hi, lo;
    tie(hi.x, lo.x) = split(dsfloat(u[0]));
    tie(hi.y, lo.y) = split(dsfloat(u[1]));
    tie(hi.z, lo.z) = split(dsfloat(u[2]));
    tie(hi.w, lo.w) = make_tuple(__int_as_float(v), 0);
    left = make_tuple(hi, lo);
}

/**
 * Pack fixed_vector<double, 3> and unsigned int into tuple of float4 and float4.
 *
 * Takes right-hand tuple of values, as returned by make_tuple.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<double, 3> >
      , boost::is_same<V, unsigned int>
    >, void>::type
operator<<=(tuple<T&, T&> left, tuple<U, V> right)
{
    fixed_vector<double, 3> u;
    unsigned int v;
    tie(u, v) = right;
    float4 hi, lo;
    tie(hi.x, lo.x) = split(dsfloat(u[0]));
    tie(hi.y, lo.y) = split(dsfloat(u[1]));
    tie(hi.z, lo.z) = split(dsfloat(u[2]));
    tie(hi.w, lo.w) = make_tuple(__int_as_float(v), 0);
    left = make_tuple(hi, lo);
}

/**
 * Unpack float4 and float4 into tuple of fixed_vector<double, 3> and unsigned int.
 *
 * Takes right-hand tuple of references, as returned by tie.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T const, float4 const>
      , boost::is_same<U, fixed_vector<double, 3> >
      , boost::is_same<V, unsigned int>
    >, void>::type
operator<<=(tuple<U&, V&> left, tuple<T&, T&> right)
{
    float4 hi, lo;
    tie(hi, lo) = right;
    fixed_vector<double, 3> u;
    u[0] = dsfloat(hi.x, lo.x);
    u[1] = dsfloat(hi.y, lo.y);
    u[2] = dsfloat(hi.z, lo.z);
    unsigned int v = __float_as_int(hi.w);
    left = make_tuple(u, v);
}

/**
 * Unpack float4 and float4 into tuple of fixed_vector<double, 3> and unsigned int.
 *
 * Takes right-hand tuple of values, as returned by make_tuple.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<double, 3> >
      , boost::is_same<V, unsigned int>
    >, void>::type
operator<<=(tuple<U&, V&> left, tuple<T, T> right)
{
    float4 hi, lo;
    tie(hi, lo) = right;
    fixed_vector<double, 3> u;
    u[0] = dsfloat(hi.x, lo.x);
    u[1] = dsfloat(hi.y, lo.y);
    u[2] = dsfloat(hi.z, lo.z);
    unsigned int v = __float_as_int(hi.w);
    left = make_tuple(u, v);
}

/**
 * Pack fixed_vector<double, 2> and unsigned int into tuple of float4 and float4.
 *
 * Takes right-hand tuple of references, as returned by tie.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U const, fixed_vector<double, 2> const>
      , boost::is_same<V const, unsigned int const>
    >, void>::type
operator<<=(tuple<T&, T&> left, tuple<U&, V&> right)
{
    fixed_vector<double, 2> u;
    unsigned int v;
    tie(u, v) = right;
    float4 hi, lo;
    tie(hi.x, lo.x) = split(dsfloat(u[0]));
    tie(hi.y, lo.y) = split(dsfloat(u[1]));
    tie(hi.z, lo.z) = make_tuple(0, 0);
    tie(hi.w, lo.w) = make_tuple(__int_as_float(v), 0);
    left = make_tuple(hi, lo);
}

/**
 * Pack fixed_vector<double, 2> and unsigned int into tuple of float4 and float4.
 *
 * Takes right-hand tuple of values, as returned by make_tuple.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<double, 2> >
      , boost::is_same<V, unsigned int>
    >, void>::type
operator<<=(tuple<T&, T&> left, tuple<U, V> right)
{
    fixed_vector<double, 2> u;
    unsigned int v;
    tie(u, v) = right;
    float4 hi, lo;
    tie(hi.x, lo.x) = split(dsfloat(u[0]));
    tie(hi.y, lo.y) = split(dsfloat(u[1]));
    tie(hi.z, lo.z) = make_tuple(0, 0);
    tie(hi.w, lo.w) = make_tuple(__int_as_float(v), 0);
    left = make_tuple(hi, lo);
}

/**
 * Unpack float4 and float4 into tuple of fixed_vector<double, 2> and unsigned int.
 *
 * Takes right-hand tuple of references, as returned by tie.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T const, float4 const>
      , boost::is_same<U, fixed_vector<double, 2> >
      , boost::is_same<V, unsigned int>
    >, void>::type
operator<<=(tuple<U&, V&> left, tuple<T&, T&> right)
{
    float4 hi, lo;
    tie(hi, lo) = right;
    fixed_vector<double, 2> u;
    u[0] = dsfloat(hi.x, lo.x);
    u[1] = dsfloat(hi.y, lo.y);
    unsigned int v = __float_as_int(hi.w);
    left = make_tuple(u, v);
}

/**
 * Unpack float4 and float4 into tuple of fixed_vector<double, 2> and unsigned int.
 *
 * Takes right-hand tuple of values, as returned by make_tuple.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<double, 2> >
      , boost::is_same<V, unsigned int>
    >, void>::type
operator<<=(tuple<U&, V&> left, tuple<T, T> right)
{
    float4 hi, lo;
    tie(hi, lo) = right;
    fixed_vector<double, 2> u;
    u[0] = dsfloat(hi.x, lo.x);
    u[1] = dsfloat(hi.y, lo.y);
    unsigned int v = __float_as_int(hi.w);
    left = make_tuple(u, v);
}

/**
 * Pack fixed_vector<double, 3> and float into tuple of float4 and float4.
 *
 * Takes right-hand tuple of references, as returned by tie.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U const, fixed_vector<double, 3> const>
      , boost::is_same<V const, float const>
    >, void>::type
operator<<=(tuple<T&, T&> left, tuple<U&, V&> right)
{
    fixed_vector<double, 3> u;
    float v;
    tie(u, v) = right;
    float4 hi, lo;
    tie(hi.x, lo.x) = split(dsfloat(u[0]));
    tie(hi.y, lo.y) = split(dsfloat(u[1]));
    tie(hi.z, lo.z) = split(dsfloat(u[2]));
    tie(hi.w, lo.w) = make_tuple(v, 0);
    left = make_tuple(hi, lo);
}

/**
 * Pack fixed_vector<double, 3> and float into tuple of float4 and float4.
 *
 * Takes right-hand tuple of values, as returned by make_tuple.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<double, 3> >
      , boost::is_same<V, float>
    >, void>::type
operator<<=(tuple<T&, T&> left, tuple<U, V> right)
{
    fixed_vector<double, 3> u;
    float v;
    tie(u, v) = right;
    float4 hi, lo;
    tie(hi.x, lo.x) = split(dsfloat(u[0]));
    tie(hi.y, lo.y) = split(dsfloat(u[1]));
    tie(hi.z, lo.z) = split(dsfloat(u[2]));
    tie(hi.w, lo.w) = make_tuple(v, 0);
    left = make_tuple(hi, lo);
}

/**
 * Unpack float4 and float4 into tuple of fixed_vector<double, 3> and float.
 *
 * Takes right-hand tuple of references, as returned by tie.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T const, float4 const>
      , boost::is_same<U, fixed_vector<double, 3> >
      , boost::is_same<V, float>
    >, void>::type
operator<<=(tuple<U&, V&> left, tuple<T&, T&> right)
{
    float4 hi, lo;
    tie(hi, lo) = right;
    fixed_vector<double, 3> u;
    u[0] = dsfloat(hi.x, lo.x);
    u[1] = dsfloat(hi.y, lo.y);
    u[2] = dsfloat(hi.z, lo.z);
    float v = hi.w;
    left = make_tuple(u, v);
}

/**
 * Unpack float4 and float4 into tuple of fixed_vector<double, 3> and float.
 *
 * Takes right-hand tuple of values, as returned by make_tuple.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<double, 3> >
      , boost::is_same<V, float>
    >, void>::type
operator<<=(tuple<U&, V&> left, tuple<T, T> right)
{
    float4 hi, lo;
    tie(hi, lo) = right;
    fixed_vector<double, 3> u;
    u[0] = dsfloat(hi.x, lo.x);
    u[1] = dsfloat(hi.y, lo.y);
    u[2] = dsfloat(hi.z, lo.z);
    float v = hi.w;
    left = make_tuple(u, v);
}

/**
 * Pack fixed_vector<double, 2> and float into tuple of float4 and float4.
 *
 * Takes right-hand tuple of references, as returned by tie.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U const, fixed_vector<double, 2> const>
      , boost::is_same<V const, float const>
    >, void>::type
operator<<=(tuple<T&, T&> left, tuple<U&, V&> right)
{
    fixed_vector<double, 2> u;
    float v;
    tie(u, v) = right;
    float4 hi, lo;
    tie(hi.x, lo.x) = split(dsfloat(u[0]));
    tie(hi.y, lo.y) = split(dsfloat(u[1]));
    tie(hi.z, lo.z) = make_tuple(0, 0);
    tie(hi.w, lo.w) = make_tuple(v, 0);
    left = make_tuple(hi, lo);
}

/**
 * Pack fixed_vector<double, 2> and float into tuple of float4 and float4.
 *
 * Takes right-hand tuple of values, as returned by make_tuple.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<double, 2> >
      , boost::is_same<V, float>
    >, void>::type
operator<<=(tuple<T&, T&> left, tuple<U, V> right)
{
    fixed_vector<double, 2> u;
    float v;
    tie(u, v) = right;
    float4 hi, lo;
    tie(hi.x, lo.x) = split(dsfloat(u[0]));
    tie(hi.y, lo.y) = split(dsfloat(u[1]));
    tie(hi.z, lo.z) = make_tuple(0, 0);
    tie(hi.w, lo.w) = make_tuple(v, 0);
    left = make_tuple(hi, lo);
}

/**
 * Unpack float4 and float4 into tuple of fixed_vector<double, 2> and float.
 *
 * Takes right-hand tuple of references, as returned by tie.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T const, float4 const>
      , boost::is_same<U, fixed_vector<double, 2> >
      , boost::is_same<V, float>
    >, void>::type
operator<<=(tuple<U&, V&> left, tuple<T&, T&> right)
{
    float4 hi, lo;
    tie(hi, lo) = right;
    fixed_vector<double, 2> u;
    u[0] = dsfloat(hi.x, lo.x);
    u[1] = dsfloat(hi.y, lo.y);
    float v = hi.w;
    left = make_tuple(u, v);
}

/**
 * Unpack float4 and float4 into tuple of fixed_vector<double, 2> and float.
 *
 * Takes right-hand tuple of values, as returned by make_tuple.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<double, 2> >
      , boost::is_same<V, float>
    >, void>::type
operator<<=(tuple<U&, V&> left, tuple<T, T> right)
{
    float4 hi, lo;
    tie(hi, lo) = right;
    fixed_vector<double, 2> u;
    u[0] = dsfloat(hi.x, lo.x);
    u[1] = dsfloat(hi.y, lo.y);
    float v = hi.w;
    left = make_tuple(u, v);
}

/**
 * Pack fixed_vector<double, 3> and double into tuple of float4 and float4.
 *
 * Takes right-hand tuple of references, as returned by tie.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U const, fixed_vector<double, 3> const>
      , boost::is_same<V const, double const>
    >, void>::type
operator<<=(tuple<T&, T&> left, tuple<U&, V&> right)
{
    fixed_vector<double, 3> u;
    double v;
    tie(u, v) = right;
    float4 hi, lo;
    tie(hi.x, lo.x) = split(dsfloat(u[0]));
    tie(hi.y, lo.y) = split(dsfloat(u[1]));
    tie(hi.z, lo.z) = split(dsfloat(u[2]));
    tie(hi.w, lo.w) = split(dsfloat(v));
    left = make_tuple(hi, lo);
}

/**
 * Pack fixed_vector<double, 3> and double into tuple of float4 and float4.
 *
 * Takes right-hand tuple of values, as returned by make_tuple.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<double, 3> >
      , boost::is_same<V, double>
    >, void>::type
operator<<=(tuple<T&, T&> left, tuple<U, V> right)
{
    fixed_vector<double, 3> u;
    double v;
    tie(u, v) = right;
    float4 hi, lo;
    tie(hi.x, lo.x) = split(dsfloat(u[0]));
    tie(hi.y, lo.y) = split(dsfloat(u[1]));
    tie(hi.z, lo.z) = split(dsfloat(u[2]));
    tie(hi.w, lo.w) = split(dsfloat(v));
    left = make_tuple(hi, lo);
}

/**
 * Unpack float4 and float4 into tuple of fixed_vector<double, 3> and double.
 *
 * Takes right-hand tuple of references, as returned by tie.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T const, float4 const>
      , boost::is_same<U, fixed_vector<double, 3> >
      , boost::is_same<V, double>
    >, void>::type
operator<<=(tuple<U&, V&> left, tuple<T&, T&> right)
{
    float4 hi, lo;
    tie(hi, lo) = right;
    fixed_vector<double, 3> u;
    u[0] = dsfloat(hi.x, lo.x);
    u[1] = dsfloat(hi.y, lo.y);
    u[2] = dsfloat(hi.z, lo.z);
    double v = dsfloat(hi.w, lo.w);
    left = make_tuple(u, v);
}

/**
 * Unpack float4 and float4 into tuple of fixed_vector<double, 3> and double.
 *
 * Takes right-hand tuple of values, as returned by make_tuple.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<double, 3> >
      , boost::is_same<V, double>
    >, void>::type
operator<<=(tuple<U&, V&> left, tuple<T, T> right)
{
    float4 hi, lo;
    tie(hi, lo) = right;
    fixed_vector<double, 3> u;
    u[0] = dsfloat(hi.x, lo.x);
    u[1] = dsfloat(hi.y, lo.y);
    u[2] = dsfloat(hi.z, lo.z);
    double v = dsfloat(hi.w, lo.w);
    left = make_tuple(u, v);
}

/**
 * Pack fixed_vector<double, 2> and double into tuple of float4 and float4.
 *
 * Takes right-hand tuple of references, as returned by tie.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U const, fixed_vector<double, 2> const>
      , boost::is_same<V const, double const>
    >, void>::type
operator<<=(tuple<T&, T&> left, tuple<U&, V&> right)
{
    fixed_vector<double, 2> u;
    double v;
    tie(u, v) = right;
    float4 hi, lo;
    tie(hi.x, lo.x) = split(dsfloat(u[0]));
    tie(hi.y, lo.y) = split(dsfloat(u[1]));
    tie(hi.z, lo.z) = make_tuple(0, 0);
    tie(hi.w, lo.w) = split(dsfloat(v));
    left = make_tuple(hi, lo);
}

/**
 * Pack fixed_vector<double, 2> and double into tuple of float4 and float4.
 *
 * Takes right-hand tuple of values, as returned by make_tuple.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<double, 2> >
      , boost::is_same<V, double>
    >, void>::type
operator<<=(tuple<T&, T&> left, tuple<U, V> right)
{
    fixed_vector<double, 2> u;
    double v;
    tie(u, v) = right;
    float4 hi, lo;
    tie(hi.x, lo.x) = split(dsfloat(u[0]));
    tie(hi.y, lo.y) = split(dsfloat(u[1]));
    tie(hi.z, lo.z) = make_tuple(0, 0);
    tie(hi.w, lo.w) = split(dsfloat(v));
    left = make_tuple(hi, lo);
}

/**
 * Unpack float4 and float4 into tuple of fixed_vector<double, 2> and double.
 *
 * Takes right-hand tuple of references, as returned by tie.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T const, float4 const>
      , boost::is_same<U, fixed_vector<double, 2> >
      , boost::is_same<V, double>
    >, void>::type
operator<<=(tuple<U&, V&> left, tuple<T&, T&> right)
{
    float4 hi, lo;
    tie(hi, lo) = right;
    fixed_vector<double, 2> u;
    u[0] = dsfloat(hi.x, lo.x);
    u[1] = dsfloat(hi.y, lo.y);
    double v = dsfloat(hi.w, lo.w);
    left = make_tuple(u, v);
}

/**
 * Unpack float4 and float4 into tuple of fixed_vector<double, 2> and double.
 *
 * Takes right-hand tuple of values, as returned by make_tuple.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<double, 2> >
      , boost::is_same<V, double>
    >, void>::type
operator<<=(tuple<U&, V&> left, tuple<T, T> right)
{
    float4 hi, lo;
    tie(hi, lo) = right;
    fixed_vector<double, 2> u;
    u[0] = dsfloat(hi.x, lo.x);
    u[1] = dsfloat(hi.y, lo.y);
    double v = dsfloat(hi.w, lo.w);
    left = make_tuple(u, v);
}

/**
 * Pack fixed_vector<dsfloat, 3> and int into tuple of float4 and float4.
 *
 * Takes right-hand tuple of references, as returned by tie.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U const, fixed_vector<dsfloat, 3> const>
      , boost::is_same<V const, int const>
    >, void>::type
operator<<=(tuple<T&, T&> left, tuple<U&, V&> right)
{
    fixed_vector<dsfloat, 3> u;
    int v;
    tie(u, v) = right;
    float4 hi, lo;
    tie(hi.x, lo.x) = split(u[0]);
    tie(hi.y, lo.y) = split(u[1]);
    tie(hi.z, lo.z) = split(u[2]);
    tie(hi.w, lo.w) = make_tuple(__int_as_float(v), 0);
    left = make_tuple(hi, lo);
}

/**
 * Pack fixed_vector<dsfloat, 3> and int into tuple of float4 and float4.
 *
 * Takes right-hand tuple of values, as returned by make_tuple.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<dsfloat, 3> >
      , boost::is_same<V, int>
    >, void>::type
operator<<=(tuple<T&, T&> left, tuple<U, V> right)
{
    fixed_vector<dsfloat, 3> u;
    int v;
    tie(u, v) = right;
    float4 hi, lo;
    tie(hi.x, lo.x) = split(u[0]);
    tie(hi.y, lo.y) = split(u[1]);
    tie(hi.z, lo.z) = split(u[2]);
    tie(hi.w, lo.w) = make_tuple(__int_as_float(v), 0);
    left = make_tuple(hi, lo);
}

/**
 * Unpack float4 and float4 into tuple of fixed_vector<dsfloat, 3> and int.
 *
 * Takes right-hand tuple of references, as returned by tie.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T const, float4 const>
      , boost::is_same<U, fixed_vector<dsfloat, 3> >
      , boost::is_same<V, int>
    >, void>::type
operator<<=(tuple<U&, V&> left, tuple<T&, T&> right)
{
    float4 hi, lo;
    tie(hi, lo) = right;
    fixed_vector<dsfloat, 3> u;
    u[0] = dsfloat(hi.x, lo.x);
    u[1] = dsfloat(hi.y, lo.y);
    u[2] = dsfloat(hi.z, lo.z);
    int v = __float_as_int(hi.w);
    left = make_tuple(u, v);
}

/**
 * Unpack float4 and float4 into tuple of fixed_vector<dsfloat, 3> and int.
 *
 * Takes right-hand tuple of values, as returned by make_tuple.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<dsfloat, 3> >
      , boost::is_same<V, int>
    >, void>::type
operator<<=(tuple<U&, V&> left, tuple<T, T> right)
{
    float4 hi, lo;
    tie(hi, lo) = right;
    fixed_vector<dsfloat, 3> u;
    u[0] = dsfloat(hi.x, lo.x);
    u[1] = dsfloat(hi.y, lo.y);
    u[2] = dsfloat(hi.z, lo.z);
    int v = __float_as_int(hi.w);
    left = make_tuple(u, v);
}

/**
 * Pack fixed_vector<dsfloat, 2> and int into tuple of float4 and float4.
 *
 * Takes right-hand tuple of references, as returned by tie.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U const, fixed_vector<dsfloat, 2> const>
      , boost::is_same<V const, int const>
    >, void>::type
operator<<=(tuple<T&, T&> left, tuple<U&, V&> right)
{
    fixed_vector<dsfloat, 2> u;
    int v;
    tie(u, v) = right;
    float4 hi, lo;
    tie(hi.x, lo.x) = split(u[0]);
    tie(hi.y, lo.y) = split(u[1]);
    tie(hi.z, lo.z) = make_tuple(0, 0);
    tie(hi.w, lo.w) = make_tuple(__int_as_float(v), 0);
    left = make_tuple(hi, lo);
}

/**
 * Pack fixed_vector<dsfloat, 2> and int into tuple of float4 and float4.
 *
 * Takes right-hand tuple of values, as returned by make_tuple.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<dsfloat, 2> >
      , boost::is_same<V, int>
    >, void>::type
operator<<=(tuple<T&, T&> left, tuple<U, V> right)
{
    fixed_vector<dsfloat, 2> u;
    int v;
    tie(u, v) = right;
    float4 hi, lo;
    tie(hi.x, lo.x) = split(u[0]);
    tie(hi.y, lo.y) = split(u[1]);
    tie(hi.z, lo.z) = make_tuple(0, 0);
    tie(hi.w, lo.w) = make_tuple(__int_as_float(v), 0);
    left = make_tuple(hi, lo);
}

/**
 * Unpack float4 and float4 into tuple of fixed_vector<dsfloat, 2> and int.
 *
 * Takes right-hand tuple of references, as returned by tie.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T const, float4 const>
      , boost::is_same<U, fixed_vector<dsfloat, 2> >
      , boost::is_same<V, int>
    >, void>::type
operator<<=(tuple<U&, V&> left, tuple<T&, T&> right)
{
    float4 hi, lo;
    tie(hi, lo) = right;
    fixed_vector<dsfloat, 2> u;
    u[0] = dsfloat(hi.x, lo.x);
    u[1] = dsfloat(hi.y, lo.y);
    int v = __float_as_int(hi.w);
    left = make_tuple(u, v);
}

/**
 * Unpack float4 and float4 into tuple of fixed_vector<dsfloat, 2> and int.
 *
 * Takes right-hand tuple of values, as returned by make_tuple.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<dsfloat, 2> >
      , boost::is_same<V, int>
    >, void>::type
operator<<=(tuple<U&, V&> left, tuple<T, T> right)
{
    float4 hi, lo;
    tie(hi, lo) = right;
    fixed_vector<dsfloat, 2> u;
    u[0] = dsfloat(hi.x, lo.x);
    u[1] = dsfloat(hi.y, lo.y);
    int v = __float_as_int(hi.w);
    left = make_tuple(u, v);
}

/**
 * Pack fixed_vector<dsfloat, 3> and unsigned int into tuple of float4 and float4.
 *
 * Takes right-hand tuple of references, as returned by tie.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U const, fixed_vector<dsfloat, 3> const>
      , boost::is_same<V const, unsigned int const>
    >, void>::type
operator<<=(tuple<T&, T&> left, tuple<U&, V&> right)
{
    fixed_vector<dsfloat, 3> u;
    unsigned int v;
    tie(u, v) = right;
    float4 hi, lo;
    tie(hi.x, lo.x) = split(u[0]);
    tie(hi.y, lo.y) = split(u[1]);
    tie(hi.z, lo.z) = split(u[2]);
    tie(hi.w, lo.w) = make_tuple(__int_as_float(v), 0);
    left = make_tuple(hi, lo);
}

/**
 * Pack fixed_vector<dsfloat, 3> and unsigned int into tuple of float4 and float4.
 *
 * Takes right-hand tuple of values, as returned by make_tuple.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<dsfloat, 3> >
      , boost::is_same<V, unsigned int>
    >, void>::type
operator<<=(tuple<T&, T&> left, tuple<U, V> right)
{
    fixed_vector<dsfloat, 3> u;
    unsigned int v;
    tie(u, v) = right;
    float4 hi, lo;
    tie(hi.x, lo.x) = split(u[0]);
    tie(hi.y, lo.y) = split(u[1]);
    tie(hi.z, lo.z) = split(u[2]);
    tie(hi.w, lo.w) = make_tuple(__int_as_float(v), 0);
    left = make_tuple(hi, lo);
}

/**
 * Unpack float4 and float4 into tuple of fixed_vector<dsfloat, 3> and unsigned int.
 *
 * Takes right-hand tuple of references, as returned by tie.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T const, float4 const>
      , boost::is_same<U, fixed_vector<dsfloat, 3> >
      , boost::is_same<V, unsigned int>
    >, void>::type
operator<<=(tuple<U&, V&> left, tuple<T&, T&> right)
{
    float4 hi, lo;
    tie(hi, lo) = right;
    fixed_vector<dsfloat, 3> u;
    u[0] = dsfloat(hi.x, lo.x);
    u[1] = dsfloat(hi.y, lo.y);
    u[2] = dsfloat(hi.z, lo.z);
    unsigned int v = __float_as_int(hi.w);
    left = make_tuple(u, v);
}

/**
 * Unpack float4 and float4 into tuple of fixed_vector<dsfloat, 3> and unsigned int.
 *
 * Takes right-hand tuple of values, as returned by make_tuple.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<dsfloat, 3> >
      , boost::is_same<V, unsigned int>
    >, void>::type
operator<<=(tuple<U&, V&> left, tuple<T, T> right)
{
    float4 hi, lo;
    tie(hi, lo) = right;
    fixed_vector<dsfloat, 3> u;
    u[0] = dsfloat(hi.x, lo.x);
    u[1] = dsfloat(hi.y, lo.y);
    u[2] = dsfloat(hi.z, lo.z);
    unsigned int v = __float_as_int(hi.w);
    left = make_tuple(u, v);
}

/**
 * Pack fixed_vector<dsfloat, 2> and unsigned int into tuple of float4 and float4.
 *
 * Takes right-hand tuple of references, as returned by tie.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U const, fixed_vector<dsfloat, 2> const>
      , boost::is_same<V const, unsigned int const>
    >, void>::type
operator<<=(tuple<T&, T&> left, tuple<U&, V&> right)
{
    fixed_vector<dsfloat, 2> u;
    unsigned int v;
    tie(u, v) = right;
    float4 hi, lo;
    tie(hi.x, lo.x) = split(u[0]);
    tie(hi.y, lo.y) = split(u[1]);
    tie(hi.z, lo.z) = make_tuple(0, 0);
    tie(hi.w, lo.w) = make_tuple(__int_as_float(v), 0);
    left = make_tuple(hi, lo);
}

/**
 * Pack fixed_vector<dsfloat, 2> and unsigned int into tuple of float4 and float4.
 *
 * Takes right-hand tuple of values, as returned by make_tuple.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<dsfloat, 2> >
      , boost::is_same<V, unsigned int>
    >, void>::type
operator<<=(tuple<T&, T&> left, tuple<U, V> right)
{
    fixed_vector<dsfloat, 2> u;
    unsigned int v;
    tie(u, v) = right;
    float4 hi, lo;
    tie(hi.x, lo.x) = split(u[0]);
    tie(hi.y, lo.y) = split(u[1]);
    tie(hi.z, lo.z) = make_tuple(0, 0);
    tie(hi.w, lo.w) = make_tuple(__int_as_float(v), 0);
    left = make_tuple(hi, lo);
}

/**
 * Unpack float4 and float4 into tuple of fixed_vector<dsfloat, 2> and unsigned int.
 *
 * Takes right-hand tuple of references, as returned by tie.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T const, float4 const>
      , boost::is_same<U, fixed_vector<dsfloat, 2> >
      , boost::is_same<V, unsigned int>
    >, void>::type
operator<<=(tuple<U&, V&> left, tuple<T&, T&> right)
{
    float4 hi, lo;
    tie(hi, lo) = right;
    fixed_vector<dsfloat, 2> u;
    u[0] = dsfloat(hi.x, lo.x);
    u[1] = dsfloat(hi.y, lo.y);
    unsigned int v = __float_as_int(hi.w);
    left = make_tuple(u, v);
}

/**
 * Unpack float4 and float4 into tuple of fixed_vector<dsfloat, 2> and unsigned int.
 *
 * Takes right-hand tuple of values, as returned by make_tuple.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<dsfloat, 2> >
      , boost::is_same<V, unsigned int>
    >, void>::type
operator<<=(tuple<U&, V&> left, tuple<T, T> right)
{
    float4 hi, lo;
    tie(hi, lo) = right;
    fixed_vector<dsfloat, 2> u;
    u[0] = dsfloat(hi.x, lo.x);
    u[1] = dsfloat(hi.y, lo.y);
    unsigned int v = __float_as_int(hi.w);
    left = make_tuple(u, v);
}

/**
 * Pack fixed_vector<dsfloat, 3> and float into tuple of float4 and float4.
 *
 * Takes right-hand tuple of references, as returned by tie.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U const, fixed_vector<dsfloat, 3> const>
      , boost::is_same<V const, float const>
    >, void>::type
operator<<=(tuple<T&, T&> left, tuple<U&, V&> right)
{
    fixed_vector<dsfloat, 3> u;
    float v;
    tie(u, v) = right;
    float4 hi, lo;
    tie(hi.x, lo.x) = split(u[0]);
    tie(hi.y, lo.y) = split(u[1]);
    tie(hi.z, lo.z) = split(u[2]);
    tie(hi.w, lo.w) = make_tuple(v, 0);
    left = make_tuple(hi, lo);
}

/**
 * Pack fixed_vector<dsfloat, 3> and float into tuple of float4 and float4.
 *
 * Takes right-hand tuple of values, as returned by make_tuple.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<dsfloat, 3> >
      , boost::is_same<V, float>
    >, void>::type
operator<<=(tuple<T&, T&> left, tuple<U, V> right)
{
    fixed_vector<dsfloat, 3> u;
    float v;
    tie(u, v) = right;
    float4 hi, lo;
    tie(hi.x, lo.x) = split(u[0]);
    tie(hi.y, lo.y) = split(u[1]);
    tie(hi.z, lo.z) = split(u[2]);
    tie(hi.w, lo.w) = make_tuple(v, 0);
    left = make_tuple(hi, lo);
}

/**
 * Unpack float4 and float4 into tuple of fixed_vector<dsfloat, 3> and float.
 *
 * Takes right-hand tuple of references, as returned by tie.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T const, float4 const>
      , boost::is_same<U, fixed_vector<dsfloat, 3> >
      , boost::is_same<V, float>
    >, void>::type
operator<<=(tuple<U&, V&> left, tuple<T&, T&> right)
{
    float4 hi, lo;
    tie(hi, lo) = right;
    fixed_vector<dsfloat, 3> u;
    u[0] = dsfloat(hi.x, lo.x);
    u[1] = dsfloat(hi.y, lo.y);
    u[2] = dsfloat(hi.z, lo.z);
    float v = hi.w;
    left = make_tuple(u, v);
}

/**
 * Unpack float4 and float4 into tuple of fixed_vector<dsfloat, 3> and float.
 *
 * Takes right-hand tuple of values, as returned by make_tuple.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<dsfloat, 3> >
      , boost::is_same<V, float>
    >, void>::type
operator<<=(tuple<U&, V&> left, tuple<T, T> right)
{
    float4 hi, lo;
    tie(hi, lo) = right;
    fixed_vector<dsfloat, 3> u;
    u[0] = dsfloat(hi.x, lo.x);
    u[1] = dsfloat(hi.y, lo.y);
    u[2] = dsfloat(hi.z, lo.z);
    float v = hi.w;
    left = make_tuple(u, v);
}

/**
 * Pack fixed_vector<dsfloat, 2> and float into tuple of float4 and float4.
 *
 * Takes right-hand tuple of references, as returned by tie.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U const, fixed_vector<dsfloat, 2> const>
      , boost::is_same<V const, float const>
    >, void>::type
operator<<=(tuple<T&, T&> left, tuple<U&, V&> right)
{
    fixed_vector<dsfloat, 2> u;
    float v;
    tie(u, v) = right;
    float4 hi, lo;
    tie(hi.x, lo.x) = split(u[0]);
    tie(hi.y, lo.y) = split(u[1]);
    tie(hi.z, lo.z) = make_tuple(0, 0);
    tie(hi.w, lo.w) = make_tuple(v, 0);
    left = make_tuple(hi, lo);
}

/**
 * Pack fixed_vector<dsfloat, 2> and float into tuple of float4 and float4.
 *
 * Takes right-hand tuple of values, as returned by make_tuple.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<dsfloat, 2> >
      , boost::is_same<V, float>
    >, void>::type
operator<<=(tuple<T&, T&> left, tuple<U, V> right)
{
    fixed_vector<dsfloat, 2> u;
    float v;
    tie(u, v) = right;
    float4 hi, lo;
    tie(hi.x, lo.x) = split(u[0]);
    tie(hi.y, lo.y) = split(u[1]);
    tie(hi.z, lo.z) = make_tuple(0, 0);
    tie(hi.w, lo.w) = make_tuple(v, 0);
    left = make_tuple(hi, lo);
}

/**
 * Unpack float4 and float4 into tuple of fixed_vector<dsfloat, 2> and float.
 *
 * Takes right-hand tuple of references, as returned by tie.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T const, float4 const>
      , boost::is_same<U, fixed_vector<dsfloat, 2> >
      , boost::is_same<V, float>
    >, void>::type
operator<<=(tuple<U&, V&> left, tuple<T&, T&> right)
{
    float4 hi, lo;
    tie(hi, lo) = right;
    fixed_vector<dsfloat, 2> u;
    u[0] = dsfloat(hi.x, lo.x);
    u[1] = dsfloat(hi.y, lo.y);
    float v = hi.w;
    left = make_tuple(u, v);
}

/**
 * Unpack float4 and float4 into tuple of fixed_vector<dsfloat, 2> and float.
 *
 * Takes right-hand tuple of values, as returned by make_tuple.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<dsfloat, 2> >
      , boost::is_same<V, float>
    >, void>::type
operator<<=(tuple<U&, V&> left, tuple<T, T> right)
{
    float4 hi, lo;
    tie(hi, lo) = right;
    fixed_vector<dsfloat, 2> u;
    u[0] = dsfloat(hi.x, lo.x);
    u[1] = dsfloat(hi.y, lo.y);
    float v = hi.w;
    left = make_tuple(u, v);
}

/**
 * Pack fixed_vector<dsfloat, 3> and dsfloat into tuple of float4 and float4.
 *
 * Takes right-hand tuple of references, as returned by tie.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U const, fixed_vector<dsfloat, 3> const>
      , boost::is_same<V const, dsfloat const>
    >, void>::type
operator<<=(tuple<T&, T&> left, tuple<U&, V&> right)
{
    fixed_vector<dsfloat, 3> u;
    dsfloat v;
    tie(u, v) = right;
    float4 hi, lo;
    tie(hi.x, lo.x) = split(u[0]);
    tie(hi.y, lo.y) = split(u[1]);
    tie(hi.z, lo.z) = split(u[2]);
    tie(hi.w, lo.w) = split(v);
    left = make_tuple(hi, lo);
}

/**
 * Pack fixed_vector<dsfloat, 3> and dsfloat into tuple of float4 and float4.
 *
 * Takes right-hand tuple of values, as returned by make_tuple.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<dsfloat, 3> >
      , boost::is_same<V, dsfloat>
    >, void>::type
operator<<=(tuple<T&, T&> left, tuple<U, V> right)
{
    fixed_vector<dsfloat, 3> u;
    dsfloat v;
    tie(u, v) = right;
    float4 hi, lo;
    tie(hi.x, lo.x) = split(u[0]);
    tie(hi.y, lo.y) = split(u[1]);
    tie(hi.z, lo.z) = split(u[2]);
    tie(hi.w, lo.w) = split(v);
    left = make_tuple(hi, lo);
}

/**
 * Unpack float4 and float4 into tuple of fixed_vector<dsfloat, 3> and dsfloat.
 *
 * Takes right-hand tuple of references, as returned by tie.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T const, float4 const>
      , boost::is_same<U, fixed_vector<dsfloat, 3> >
      , boost::is_same<V, dsfloat>
    >, void>::type
operator<<=(tuple<U&, V&> left, tuple<T&, T&> right)
{
    float4 hi, lo;
    tie(hi, lo) = right;
    fixed_vector<dsfloat, 3> u;
    u[0] = dsfloat(hi.x, lo.x);
    u[1] = dsfloat(hi.y, lo.y);
    u[2] = dsfloat(hi.z, lo.z);
    dsfloat v = dsfloat(hi.w, lo.w);
    left = make_tuple(u, v);
}

/**
 * Unpack float4 and float4 into tuple of fixed_vector<dsfloat, 3> and dsfloat.
 *
 * Takes right-hand tuple of values, as returned by make_tuple.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<dsfloat, 3> >
      , boost::is_same<V, dsfloat>
    >, void>::type
operator<<=(tuple<U&, V&> left, tuple<T, T> right)
{
    float4 hi, lo;
    tie(hi, lo) = right;
    fixed_vector<dsfloat, 3> u;
    u[0] = dsfloat(hi.x, lo.x);
    u[1] = dsfloat(hi.y, lo.y);
    u[2] = dsfloat(hi.z, lo.z);
    dsfloat v = dsfloat(hi.w, lo.w);
    left = make_tuple(u, v);
}

/**
 * Pack fixed_vector<dsfloat, 2> and dsfloat into tuple of float4 and float4.
 *
 * Takes right-hand tuple of references, as returned by tie.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U const, fixed_vector<dsfloat, 2> const>
      , boost::is_same<V const, dsfloat const>
    >, void>::type
operator<<=(tuple<T&, T&> left, tuple<U&, V&> right)
{
    fixed_vector<dsfloat, 2> u;
    dsfloat v;
    tie(u, v) = right;
    float4 hi, lo;
    tie(hi.x, lo.x) = split(u[0]);
    tie(hi.y, lo.y) = split(u[1]);
    tie(hi.z, lo.z) = make_tuple(0, 0);
    tie(hi.w, lo.w) = split(v);
    left = make_tuple(hi, lo);
}

/**
 * Pack fixed_vector<dsfloat, 2> and dsfloat into tuple of float4 and float4.
 *
 * Takes right-hand tuple of values, as returned by make_tuple.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<dsfloat, 2> >
      , boost::is_same<V, dsfloat>
    >, void>::type
operator<<=(tuple<T&, T&> left, tuple<U, V> right)
{
    fixed_vector<dsfloat, 2> u;
    dsfloat v;
    tie(u, v) = right;
    float4 hi, lo;
    tie(hi.x, lo.x) = split(u[0]);
    tie(hi.y, lo.y) = split(u[1]);
    tie(hi.z, lo.z) = make_tuple(0, 0);
    tie(hi.w, lo.w) = split(v);
    left = make_tuple(hi, lo);
}

/**
 * Unpack float4 and float4 into tuple of fixed_vector<dsfloat, 2> and dsfloat.
 *
 * Takes right-hand tuple of references, as returned by tie.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T const, float4 const>
      , boost::is_same<U, fixed_vector<dsfloat, 2> >
      , boost::is_same<V, dsfloat>
    >, void>::type
operator<<=(tuple<U&, V&> left, tuple<T&, T&> right)
{
    float4 hi, lo;
    tie(hi, lo) = right;
    fixed_vector<dsfloat, 2> u;
    u[0] = dsfloat(hi.x, lo.x);
    u[1] = dsfloat(hi.y, lo.y);
    dsfloat v = dsfloat(hi.w, lo.w);
    left = make_tuple(u, v);
}

/**
 * Unpack float4 and float4 into tuple of fixed_vector<dsfloat, 2> and dsfloat.
 *
 * Takes right-hand tuple of values, as returned by make_tuple.
 */
template <typename T, typename U, typename V>
inline HALMD_GPU_ENABLED typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<T, float4>
      , boost::is_same<U, fixed_vector<dsfloat, 2> >
      , boost::is_same<V, dsfloat>
    >, void>::type
operator<<=(tuple<U&, V&> left, tuple<T, T> right)
{
    float4 hi, lo;
    tie(hi, lo) = right;
    fixed_vector<dsfloat, 2> u;
    u[0] = dsfloat(hi.x, lo.x);
    u[1] = dsfloat(hi.y, lo.y);
    dsfloat v = dsfloat(hi.w, lo.w);
    left = make_tuple(u, v);
}

} // namespace detail
} // namespace blas
} // namespace numeric
} // namespace halmd

#endif /* ! HALMD_NUMERIC_BLAS_DETAIL_CUDA_VECTOR_CONVERTER_HPP */
