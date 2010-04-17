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

#include <boost/utility/enable_if.hpp>

namespace halmd { namespace mdsim { namespace gpu { namespace particle_kernel
{

/** placeholder particle */
enum { PLACEHOLDER = -1U };

/**
 * Convert particle position and tag to coalesced vector type
 */
template <typename vector_type>
__device__ inline typename boost::enable_if_c<vector_type::static_size == 3, float4>::type
tagged(vector_type v, unsigned int tag)
{
    return make_float4(v[0], v[1], v[2], __int_as_float(tag));
}

template <typename vector_type>
__device__ inline typename boost::enable_if_c<vector_type::static_size == 2, float4>::type
tagged(vector_type v, unsigned int tag)
{
    return make_float4(v[0], v[1], 0, __int_as_float(tag));
}

/**
 * Convert coalesced vector type to particle position and tag
 */
template <typename vector_type>
__device__ inline typename boost::enable_if_c<vector_type::static_size == 3, vector_type>::type
untagged(float4 v, unsigned int& tag)
{
    vector_type w;
    w[0] = v.x;
    w[1] = v.y;
    w[2] = v.z;
    tag = __float_as_int(v.w);
    return w;
}

template <typename vector_type>
__device__ inline typename boost::enable_if_c<vector_type::static_size == 2, vector_type>::type
untagged(float4 v, unsigned int& tag)
{
    vector_type w;
    w[0] = v.x;
    w[1] = v.y;
    tag = __float_as_int(v.w);
    return w;
}

/**
 * Convert coalesced vector type to particle position
 */
template <typename vector_type>
__device__ inline typename boost::enable_if_c<vector_type::static_size == 3, vector_type>::type
untagged(float4 v)
{
    vector_type w;
    w[0] = v.x;
    w[1] = v.y;
    w[2] = v.z;
    return w;
}

template <typename vector_type>
__device__ inline typename boost::enable_if_c<vector_type::static_size == 2, vector_type>::type
untagged(float4 v)
{
    vector_type w;
    w[0] = v.x;
    w[1] = v.y;
    return w;
}

}}}} // namespace halmd::mdsim::gpu::particle_kernel

#endif /* ! HALMD_MDSIM_GPU_PARTICLE_KERNEL_CUH */
