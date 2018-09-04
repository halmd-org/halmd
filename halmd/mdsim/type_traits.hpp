/*
 * Copyright © 2010  Felix Höfling
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

#ifndef HALMD_MDSIM_TYPE_TRAITS_HPP
#define HALMD_MDSIM_TYPE_TRAITS_HPP

#include <halmd/config.hpp>

#ifdef HALMD_WITH_GPU
#include <cuda_wrapper/cuda_wrapper.hpp>
#endif

#include <halmd/numeric/blas/blas.hpp>

namespace halmd {
namespace mdsim {

#ifdef HALMD_WITH_GPU
namespace detail {
namespace gpu {

/**
 * basic GPU vector types
 * e.g., float implies float2, float3, float4
 */

template <int N, typename value_type>
struct basic_vector_type;

template <> struct basic_vector_type<1, float> { typedef float type; typedef type* ptr_type; typedef type const* const_ptr_type; };
template <> struct basic_vector_type<2, float> { typedef float2 type; typedef type* ptr_type; typedef type const* const_ptr_type; };
template <> struct basic_vector_type<3, float> { typedef float3 type; typedef type* ptr_type; typedef type const* const_ptr_type; };
template <> struct basic_vector_type<4, float> { typedef float4 type; typedef type* ptr_type; typedef type const* const_ptr_type; };
template <> struct basic_vector_type<6, float> { typedef fixed_vector<float, 6> type; typedef type* ptr_type; typedef type const* const_ptr_type; };

template <> struct basic_vector_type<1, dsfloat> { typedef dsfloat type; typedef dsfloat_ptr<float> ptr_type; typedef dsfloat_const_ptr<float> const_ptr_type; };
template <> struct basic_vector_type<2, dsfloat> { typedef fixed_vector<dsfloat, 2> type; typedef dsfloat_ptr<float2> ptr_type; typedef dsfloat_const_ptr<float2> const_ptr_type; };
template <> struct basic_vector_type<3, dsfloat> { typedef fixed_vector<dsfloat, 3> type; typedef dsfloat_ptr<float3> ptr_type; typedef dsfloat_const_ptr<float3> const_ptr_type; };
template <> struct basic_vector_type<4, dsfloat> { typedef fixed_vector<dsfloat, 4> type; typedef dsfloat_ptr<float4> ptr_type; typedef dsfloat_const_ptr<float4> const_ptr_type; };
template <> struct basic_vector_type<6, dsfloat> { typedef fixed_vector<dsfloat, 6> type; typedef dsfloat_ptr<type> ptr_type; typedef dsfloat_const_ptr<type> const_ptr_type; };

// template <> struct basic_vector_type<1, double> { typedef double type; typedef type const* const_ptr_type; };
// template <> struct basic_vector_type<2, double> { typedef double2 type; typedef type const* const_ptr_type; };
// template <> struct basic_vector_type<3, double> { typedef double3 type; typedef type const* const_ptr_type; };
// template <> struct basic_vector_type<4, double> { typedef double4 type; typedef type const* const_ptr_type; };

template <> struct basic_vector_type<1, int> { typedef int type; typedef type* ptr_type; typedef type const* const_ptr_type; };
template <> struct basic_vector_type<2, int> { typedef int2 type; typedef type* ptr_type; typedef type const* const_ptr_type; };
template <> struct basic_vector_type<3, int> { typedef int3 type; typedef type* ptr_type; typedef type const* const_ptr_type; };
template <> struct basic_vector_type<4, int> { typedef int4 type; typedef type* ptr_type; typedef type const* const_ptr_type; };

template <> struct basic_vector_type<1, unsigned int> { typedef unsigned int type; typedef type* ptr_type; typedef type const* const_ptr_type; };
template <> struct basic_vector_type<2, unsigned int> { typedef uint2 type; typedef type* ptr_type; typedef type const* const_ptr_type; };
template <> struct basic_vector_type<3, unsigned int> { typedef uint3 type; typedef type* ptr_type; typedef type const* const_ptr_type; };
template <> struct basic_vector_type<4, unsigned int> { typedef uint4 type; typedef type* ptr_type; typedef type const* const_ptr_type; };

/**
 * definition of non-scalar GPU types
 */

template <int dimension, typename value_type>
struct type_traits;

template <typename value_type>
struct type_traits<4, value_type>
{
    typedef typename basic_vector_type<4, value_type>::type vector_type;
    typedef typename basic_vector_type<4, value_type>::type coalesced_vector_type;
    typedef typename basic_vector_type<6, value_type>::type coalesced_pseudo_vector_type;
    typedef typename basic_vector_type<4, value_type>::ptr_type ptr_type;
    typedef typename basic_vector_type<4, value_type>::const_ptr_type const_ptr_type;
};

template <typename value_type>
struct type_traits<3, value_type>
{
    typedef typename basic_vector_type<3, value_type>::type vector_type;
    typedef typename basic_vector_type<4, value_type>::type coalesced_vector_type;
    typedef typename basic_vector_type<4, value_type>::type coalesced_pseudo_vector_type;
    typedef typename basic_vector_type<4, value_type>::ptr_type ptr_type;
    typedef typename basic_vector_type<4, value_type>::const_ptr_type const_ptr_type;
};

template <typename value_type>
struct type_traits<2, value_type>
{
    typedef typename basic_vector_type<2, value_type>::type vector_type;
    typedef typename basic_vector_type<2, value_type>::type coalesced_vector_type;
    typedef value_type coalesced_pseudo_vector_type;
    typedef typename basic_vector_type<2, value_type>::ptr_type ptr_type;
    typedef typename basic_vector_type<2, value_type>::const_ptr_type const_ptr_type;
};

template <typename value_type>
struct type_traits<1, value_type>
{
    typedef value_type vector_type;
    typedef value_type coalesced_vector_type;
    typedef value_type coalesced_pseudo_vector_type;
    typedef typename basic_vector_type<1, value_type>::ptr_type ptr_type;
    typedef typename basic_vector_type<1, value_type>::const_ptr_type const_ptr_type;
};

}} // namespace detail::gpu

#endif // HALMD_WITH_GPU

/**
 * definitions of non-scalar types used by the simulation
 *
 * template parameters are space dimension and value type
 */

// struct specialization is necessary because there is no analytical formula for pseudovetor dimensionality
template <int dimension, typename value_type>
struct type_traits;

template <typename value_type>
struct type_traits<1, value_type>
{
    typedef value_type vector_type;
    typedef value_type stress_tensor_type;

#ifdef HALMD_WITH_GPU
    typedef detail::gpu::type_traits<1, value_type> gpu;
#endif
};

template <typename value_type>
struct type_traits<2, value_type>
{
    typedef fixed_vector<value_type, 2> vector_type;
    // diagonal elements are the first 'dimension' elements
    // followed by the off-diagonals T_12, T_13, ..., T_23, ...
    typedef fixed_vector<value_type, 3> stress_tensor_type;
    typedef value_type pseudo_vector_type;

#ifdef HALMD_WITH_GPU
    typedef detail::gpu::type_traits<2, value_type> gpu;
#endif
};

template <typename value_type>
struct type_traits<3, value_type>
{
    typedef fixed_vector<value_type, 3> vector_type;
    // diagonal elements are the first 'dimension' elements
    // followed by the off-diagonals T_12, T_13, ..., T_23, ...
    typedef fixed_vector<value_type, 6> stress_tensor_type;
    typedef fixed_vector<value_type, 3> pseudo_vector_type;


#ifdef HALMD_WITH_GPU
    typedef detail::gpu::type_traits<3, value_type> gpu;
#endif
};

template <typename value_type>
struct type_traits<4, value_type>
{
    typedef fixed_vector<value_type, 4> vector_type;
    // diagonal elements are the first 'dimension' elements
    // followed by the off-diagonals T_12, T_13, ..., T_23, ...
    typedef fixed_vector<value_type, 10> stress_tensor_type;
    typedef fixed_vector<value_type, 6> pseudo_vector_type;

#ifdef HALMD_WITH_GPU
    typedef detail::gpu::type_traits<4, value_type> gpu;
#endif
};

} // namespace mdsim
} // namespace halmd

#endif /* ! HALMD_MDSIM_TYPE_TRAITS_HPP */
