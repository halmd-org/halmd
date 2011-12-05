/*
 * Copyright © 2008-2011  Peter Colberg, Felix Höfling, and Michael Kopp
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

namespace halmd {
namespace algorithm {
namespace gpu {

using boost::disable_if;
using boost::enable_if;
using boost::is_same;

/**
 * transformation types
 */
struct identity_;
struct negate_;
struct square_;
struct sqrt_;
// FIXME template <int index> struct at_;
struct at_0;
template <int dimension>
struct trace;
struct sum_;
struct complex_sum_;
struct quaternion_sum_;
struct max_;
struct accumulate_;

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
__device__ typename enable_if<is_same<transform_, negate_>, output_type>::type
transform(input_type v)
{
    return -v;
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

template <typename transform_, typename input_type, typename output_type>
__device__ typename enable_if<is_same<transform_, at_0>, output_type>::type
transform(input_type v)
{
    return v[0];
}

/** return trace of tensor, i.e. the sum of the first few elements */
template <typename transform_, typename input_type, typename output_type>
__device__ typename enable_if<is_same<transform_, trace<2> >, output_type>::type
transform(input_type v)
{
    return v[0] + v[1];
}

template <typename transform_, typename input_type, typename output_type>
__device__ typename enable_if<is_same<transform_, trace<3> >, output_type>::type
transform(input_type v)
{
    return v[0] + v[1] + v[2];
}

template <typename transform_, typename T>
__device__ typename boost::enable_if<boost::is_same<transform_, accumulate_>, void>::type
transform(T& acc, T&, T const& acc_, T const&)
{
    acc(acc_);
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

template <typename transform_, typename T0, typename T1>
__device__ typename enable_if<is_same<transform_, complex_sum_>, void>::type
transform(T0& r1, T1& i1, T0 r2, T1 i2)
{
    r1 += r2;
    i1 += i2;
}

template <typename transform_, typename T0, typename T1, typename T2, typename T3>
__device__ typename enable_if<is_same<transform_, quaternion_sum_>, void>::type
transform(T0& r1, T1& i1, T2& j1, T3& k1, T0 r2, T1 i2, T2 j2, T3 k2)
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
    return max(v1, v2);
}

} // namespace algorithm
} // namespace gpu
} // namespace halmd

#endif /* ! HALMD_ALGORITHM_GPU_TRANSFORM_CUH */
