/* CUDA vector types
 *
 * Copyright (C) 2008  Peter Colberg
 *
 * This program is free software: you can redistribute it and/or modify
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

#ifndef CUDA_VECTOR_TYPES_HPP
#define CUDA_VECTOR_TYPES_HPP

#include <cuda/cuda_runtime.h>


namespace cuda
{

/*
 * CUDA vector types
 */
template <unsigned dimension, typename T>
struct vector_type;

template <>
struct vector_type<1, char>
{
    typedef char1 value_type;
};

template <>
struct vector_type<1, unsigned char>
{
    typedef uchar1 value_type;
};

template <>
struct vector_type<2, char>
{
    typedef char2 value_type;
};

template <>
struct vector_type<2, unsigned char>
{
    typedef uchar2 value_type;
};

template <>
struct vector_type<3, char>
{
    typedef char3 value_type;
};

template <>
struct vector_type<3, unsigned char>
{
    typedef uchar3 value_type;
};

template <>
struct vector_type<4, char>
{
    typedef char4 value_type;
};

template <>
struct vector_type<4, unsigned char>
{
    typedef uchar4 value_type;
};

template <>
struct vector_type<1, short>
{
    typedef short1 value_type;
};

template <>
struct vector_type<1, unsigned short>
{
    typedef ushort1 value_type;
};

template <>
struct vector_type<2, short>
{
    typedef short2 value_type;
};

template <>
struct vector_type<2, unsigned short>
{
    typedef ushort2 value_type;
};

template <>
struct vector_type<3, short>
{
    typedef short3 value_type;
};

template <>
struct vector_type<3, unsigned short>
{
    typedef ushort3 value_type;
};

template <>
struct vector_type<4, short>
{
    typedef short4 value_type;
};

template <>
struct vector_type<4, unsigned short>
{
    typedef ushort4 value_type;
};

template <>
struct vector_type<1, int>
{
    typedef int1 value_type;
};

template <>
struct vector_type<1, unsigned int>
{
    typedef uint1 value_type;
};

template <>
struct vector_type<2, int>
{
    typedef int2 value_type;
};

template <>
struct vector_type<2, unsigned int>
{
    typedef uint2 value_type;
};

template <>
struct vector_type<3, int>
{
    typedef int3 value_type;
};

template <>
struct vector_type<3, unsigned int>
{
    typedef uint3 value_type;
};

template <>
struct vector_type<4, int>
{
    typedef int4 value_type;
};

template <>
struct vector_type<4, unsigned int>
{
    typedef uint4 value_type;
};

template <>
struct vector_type<1, long>
{
    typedef long1 value_type;
};

template <>
struct vector_type<1, unsigned long>
{
    typedef ulong1 value_type;
};

template <>
struct vector_type<2, long>
{
    typedef long2 value_type;
};

template <>
struct vector_type<2, unsigned long>
{
    typedef ulong2 value_type;
};

template <>
struct vector_type<3, long>
{
    typedef long3 value_type;
};

template <>
struct vector_type<3, unsigned long>
{
    typedef ulong3 value_type;
};

template <>
struct vector_type<4, long>
{
    typedef long4 value_type;
};

template <>
struct vector_type<4, unsigned long>
{
    typedef ulong4 value_type;
};

template <>
struct vector_type<1, float>
{
    typedef float1 value_type;
};

template <>
struct vector_type<2, float>
{
    typedef float2 value_type;
};

template <>
struct vector_type<3, float>
{
    typedef float3 value_type;
};

template <>
struct vector_type<4, float>
{
    typedef float4 value_type;
};

} // namespace cuda

#endif /* ! CUDA_VECTOR_TYPES_HPP */
