/*
 * Copyright © 2008-2009  Peter Colberg
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

#ifndef HALMD_ALGORITHM_GPU_TRANSFORM_CUH
#define HALMD_ALGORITHM_GPU_TRANSFORM_CUH

#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>

namespace halmd
{
namespace algorithm { namespace gpu
{

using boost::disable_if;
using boost::enable_if;
using boost::is_same;

/**
 * transformation types
 */
struct identity_;
struct square_;
struct sqrt_;
struct sum_;
struct complex_sum_;
struct quaternion_sum_;
struct max_;

/**
 * unary transformations
 */
template <typename transform_, typename input_type, typename output_type>
__device__ typename enable_if<is_same<transform_, identity_>, output_type>::type
transform(input_type v)
{
    return v;
}

template <typename transform_, typename input_type, typename output_type>
__device__ typename enable_if<is_same<transform_, square_>, output_type>::type
transform(input_type v)
{
    return inner_prod(v, v);
}

template <typename transform_, typename input_type, typename output_type>
__device__ typename enable_if<is_same<transform_, sqrt_>, output_type>::type
transform(input_type v)
{
    return sqrtf(v);
}

/**
 * binary transformations
 */
template <typename transform_, typename T>
__device__ typename enable_if<is_same<transform_, sum_>, T>::type
transform(T v1, T v2)
{
    return v1 + v2;
}

template <typename transform_, typename T>
__device__ typename enable_if<is_same<transform_, complex_sum_>, void>::type
transform(T& r1, T& i1, T r2, T i2)
{
    r1 += r2;
    i1 += i2;
}

template <typename transform_, typename T>
__device__ typename enable_if<is_same<transform_, quaternion_sum_>, void>::type
transform(T& r1, T& i1, T& j1, T& k1, T r2, T i2, T j2, T k2)
{
    r1 += r2;
    i1 += i2;
    j1 += j2;
    k1 += k2;
}

template <typename transform_, typename T>
__device__ typename enable_if<is_same<transform_, max_>, T>::type
transform(T v1, T v2)
{
    return fmaxf(v1, v2);
}

}} // namespace algorithm::gpu

} // namespace halmd

#endif /* ! HALMD_ALGORITHM_GPU_TRANSFORM_CUH */
