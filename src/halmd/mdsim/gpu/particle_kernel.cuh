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

#ifndef HALMD_MDSIM_GPU_PARTICLE_KERNEL_CUH
#define HALMD_MDSIM_GPU_PARTICLE_KERNEL_CUH

#include <boost/mpl/and.hpp>
#include <boost/mpl/int.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>
#include <cuda_runtime.h> // if compiled with C++ compiler
#include <string.h>

#ifdef __CUDACC__
# include <halmd/algorithm/gpu/tuple.cuh>
#else
# include <boost/tuple/tuple.hpp>
#endif
#include <halmd/numeric/gpu/blas/dsfloat.cuh>
#include <halmd/numeric/gpu/blas/vector.cuh>

namespace halmd
{
namespace mdsim { namespace gpu { namespace particle_kernel
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
using numeric::gpu::blas::dsfloat;

/** placeholder particle */
enum { PLACEHOLDER = -1U };

/**
 * Convert particle position and tag to coalesced vector type
 */
template <typename vector_type>
__device__ typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<boost::mpl::int_<vector_type::static_size>, boost::mpl::int_<3> >
      , boost::is_same<typename vector_type::value_type, float>
    >
  , float4
>::type
tagged(vector_type v, unsigned int tag)
{
    float4 w;
    w.x = v[0];
    w.y = v[1];
    w.z = v[2];
#ifdef __CUDACC__
    w.w = __int_as_float(tag);
#else
    // We use memcpy instead of pointer casting to avoid dereferencing a
    // type-punned pointer, which would break strict-aliasing rules.
    // Same endianness is assumed on GPU and host.
    memcpy(&w.w, &tag, sizeof(tag));
#endif
    return w;
}

template <typename vector_type>
__device__ typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<boost::mpl::int_<vector_type::static_size>, boost::mpl::int_<3> >
      , boost::is_same<typename vector_type::value_type, dsfloat>
    >
  , tuple<float4, float4>
>::type
tagged(vector_type v, unsigned int tag)
{
    float4 hi, lo;
    hi.x = dsfloat_hi(v[0]);
    lo.x = dsfloat_lo(v[0]);
    hi.y = dsfloat_hi(v[1]);
    lo.y = dsfloat_lo(v[1]);
    hi.z = dsfloat_hi(v[2]);
    lo.z = dsfloat_lo(v[2]);
#ifdef __CUDACC__
    hi.w = __int_as_float(tag);
#else
    // We use memcpy instead of pointer casting to avoid dereferencing a
    // type-punned pointer, which would break strict-aliasing rules.
    // Same endianness is assumed on GPU and host.
    memcpy(&hi.w, &tag, sizeof(tag));
#endif
    lo.w = 0;
    return make_tuple(hi, lo);
}

template <typename vector_type>
__device__ typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<boost::mpl::int_<vector_type::static_size>, boost::mpl::int_<2> >
      , boost::is_same<typename vector_type::value_type, float>
    >
  , float4
>::type
tagged(vector_type v, unsigned int tag)
{
    float4 w;
    w.x = v[0];
    w.y = v[1];
    w.z = 0;
#ifdef __CUDACC__
    w.w = __int_as_float(tag);
#else
    // We use memcpy instead of pointer casting to avoid dereferencing a
    // type-punned pointer, which would break strict-aliasing rules.
    // Same endianness is assumed on GPU and host.
    memcpy(&w.w, &tag, sizeof(tag));
#endif
    return w;
}

template <typename vector_type>
__device__ typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<boost::mpl::int_<vector_type::static_size>, boost::mpl::int_<2> >
      , boost::is_same<typename vector_type::value_type, dsfloat>
    >
  , tuple<float4, float4>
>::type
tagged(vector_type v, unsigned int tag)
{
    float4 hi, lo;
    hi.x = dsfloat_hi(v[0]);
    lo.x = dsfloat_lo(v[0]);
    hi.y = dsfloat_hi(v[1]);
    lo.y = dsfloat_lo(v[1]);
    hi.z = 0;
    lo.z = 0;
#ifdef __CUDACC__
    hi.w = __int_as_float(tag);
#else
    // We use memcpy instead of pointer casting to avoid dereferencing a
    // type-punned pointer, which would break strict-aliasing rules.
    // Same endianness is assumed on GPU and host.
    memcpy(&hi.w, &tag, sizeof(tag));
#endif
    lo.w = 0;
    return make_tuple(hi, lo);
}

/**
 * Convert coalesced vector type to particle position and tag
 */
template <typename vector_type>
__device__ typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<boost::mpl::int_<vector_type::static_size>, boost::mpl::int_<3> >
      , boost::is_same<typename vector_type::value_type, float>
    >
  , tuple<vector_type, unsigned int>
>::type
untagged(float4 v)
{
    vector_type w;
    w[0] = v.x;
    w[1] = v.y;
    w[2] = v.z;
    unsigned int tag;
#ifdef __CUDACC__
    tag = __float_as_int(v.w);
#else
    // We use memcpy instead of pointer casting to avoid dereferencing a
    // type-punned pointer, which would break strict-aliasing rules.
    // Same endianness is assumed on GPU and host.
    memcpy(&tag, &v.w, sizeof(tag));
#endif
    return make_tuple(w, tag);
}

template <typename vector_type>
__device__ typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<boost::mpl::int_<vector_type::static_size>, boost::mpl::int_<3> >
      , boost::is_same<typename vector_type::value_type, dsfloat>
    >
  , tuple<vector_type, unsigned int>
>::type
untagged(float4 hi, float4 lo)
{
    typedef typename vector_type::value_type value_type;
    vector_type w;
    w[0] = value_type(hi.x, lo.x);
    w[1] = value_type(hi.y, lo.y);
    w[2] = value_type(hi.z, lo.z);
    unsigned int tag;
#ifdef __CUDACC__
    tag = __float_as_int(hi.w);
#else
    // We use memcpy instead of pointer casting to avoid dereferencing a
    // type-punned pointer, which would break strict-aliasing rules.
    // Same endianness is assumed on GPU and host.
    memcpy(&tag, &hi.w, sizeof(tag));
#endif
    return make_tuple(w, tag);
}

template <typename vector_type>
__device__ typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<boost::mpl::int_<vector_type::static_size>, boost::mpl::int_<2> >
      , boost::is_same<typename vector_type::value_type, float>
    >
  , tuple<vector_type, unsigned int>
>::type
untagged(float4 const& v)
{
    vector_type w;
    w[0] = v.x;
    w[1] = v.y;
    unsigned int tag;
#ifdef __CUDACC__
    tag = __float_as_int(v.w);
#else
    // We use memcpy instead of pointer casting to avoid dereferencing a
    // type-punned pointer, which would break strict-aliasing rules.
    // Same endianness is assumed on GPU and host.
    memcpy(&tag, &v.w, sizeof(tag));
#endif
    return make_tuple(w, tag);
}

template <typename vector_type>
__device__ typename boost::enable_if<
    boost::mpl::and_<
        boost::is_same<boost::mpl::int_<vector_type::static_size>, boost::mpl::int_<2> >
      , boost::is_same<typename vector_type::value_type, dsfloat>
    >
  , tuple<vector_type, unsigned int>
>::type
untagged(float4 const& hi, float4 const& lo)
{
    typedef typename vector_type::value_type value_type;
    vector_type w;
    w[0] = value_type(hi.x, lo.x);
    w[1] = value_type(hi.y, lo.y);
    unsigned int tag;
#ifdef __CUDACC__
    tag = __float_as_int(hi.w);
#else
    // We use memcpy instead of pointer casting to avoid dereferencing a
    // type-punned pointer, which would break strict-aliasing rules.
    // Same endianness is assumed on GPU and host.
    memcpy(&tag, &hi.w, sizeof(tag));
#endif
    return make_tuple(w, tag);
}

}}} // namespace mdsim::gpu::particle_kernel

} // namespace halmd

#endif /* ! HALMD_MDSIM_GPU_PARTICLE_KERNEL_CUH */
