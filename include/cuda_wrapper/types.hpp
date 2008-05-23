/* CUDA vector type for translation from native types
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

#ifndef CUDA_TYPES_HPP
#define CUDA_TYPES_HPP

#include <cuda/cuda_runtime.h>


namespace cuda { namespace types
{

/*
 * CUDA vector types for translation from native types
 */
template <unsigned dimension, typename T>
struct vector;

template <>
struct vector<1, char>
{
    typedef char1 type;
};

template <>
struct vector<1, unsigned char>
{
    typedef uchar1 type;
};

template <>
struct vector<2, char>
{
    typedef char2 type;
};

template <>
struct vector<2, unsigned char>
{
    typedef uchar2 type;
};

template <>
struct vector<3, char>
{
    typedef char3 type;
};

template <>
struct vector<3, unsigned char>
{
    typedef uchar3 type;
};

template <>
struct vector<4, char>
{
    typedef char4 type;
};

template <>
struct vector<4, unsigned char>
{
    typedef uchar4 type;
};

template <>
struct vector<1, short>
{
    typedef short1 type;
};

template <>
struct vector<1, unsigned short>
{
    typedef ushort1 type;
};

template <>
struct vector<2, short>
{
    typedef short2 type;
};

template <>
struct vector<2, unsigned short>
{
    typedef ushort2 type;
};

template <>
struct vector<3, short>
{
    typedef short3 type;
};

template <>
struct vector<3, unsigned short>
{
    typedef ushort3 type;
};

template <>
struct vector<4, short>
{
    typedef short4 type;
};

template <>
struct vector<4, unsigned short>
{
    typedef ushort4 type;
};

template <>
struct vector<1, int>
{
    typedef int1 type;
};

template <>
struct vector<1, unsigned int>
{
    typedef uint1 type;
};

template <>
struct vector<2, int>
{
    typedef int2 type;
};

template <>
struct vector<2, unsigned int>
{
    typedef uint2 type;
};

template <>
struct vector<3, int>
{
    typedef int3 type;
};

template <>
struct vector<3, unsigned int>
{
    typedef uint3 type;
};

template <>
struct vector<4, int>
{
    typedef int4 type;
};

template <>
struct vector<4, unsigned int>
{
    typedef uint4 type;
};

template <>
struct vector<1, long>
{
    typedef long1 type;
};

template <>
struct vector<1, unsigned long>
{
    typedef ulong1 type;
};

template <>
struct vector<2, long>
{
    typedef long2 type;
};

template <>
struct vector<2, unsigned long>
{
    typedef ulong2 type;
};

template <>
struct vector<3, long>
{
    typedef long3 type;
};

template <>
struct vector<3, unsigned long>
{
    typedef ulong3 type;
};

template <>
struct vector<4, long>
{
    typedef long4 type;
};

template <>
struct vector<4, unsigned long>
{
    typedef ulong4 type;
};

template <>
struct vector<1, float>
{
    typedef float1 type;
};

template <>
struct vector<2, float>
{
    typedef float2 type;
};

template <>
struct vector<3, float>
{
    typedef float3 type;
};

template <>
struct vector<4, float>
{
    typedef float4 type;
};

}} // namespace cuda::types

#endif /* ! CUDA_TYPES_HPP */
