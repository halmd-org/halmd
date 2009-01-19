/* GPU vector type assignment
 *
 * Copyright © 2008-2009  Peter Colberg
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

#ifndef LJGPU_MATH_GPU_PACK_CUH
#define LJGPU_MATH_GPU_PACK_CUH

#include <cuda/cuda_runtime.h>
#include <ljgpu/math/gpu/vector2d.cuh>
#include <ljgpu/math/gpu/vector3d.cuh>

namespace ljgpu { namespace cu
{

struct __ref_float2_int
{
    float &x, &y;
    int &z;

    __device__ inline __ref_float2_int(float& x, float& y, int& z)
	: x(x), y(y), z(z) {}

    __device__ inline __ref_float2_int& operator=(float3 const& v)
    {
	x = v.x;
	y = v.y;
	z = __float_as_int(v.z);
	return *this;
    }

    __device__ inline __ref_float2_int& operator=(float4 const& v)
    {
	x = v.x;
	y = v.y;
	z = __float_as_int(v.z);
	return *this;
    }

    __device__ inline operator float4() const
    {
	return make_float4(x, y, __int_as_float(z), 0);
    }

    __device__ inline operator float3() const
    {
	return make_float3(x, y, __int_as_float(z));
    }
};

struct __cref_float2_int
{
    float const &x, &y;
    int const &z;

    __device__ inline __cref_float2_int(float const& x, float const& y, int const& z)
	: x(x), y(y), z(z) {}

    __device__ inline operator float4() const
    {
	return make_float4(x, y, __int_as_float(z), 0);
    }

    __device__ inline operator float3() const
    {
	return make_float3(x, y, __int_as_float(z));
    }
};

struct __ref_float3
{
    float &x, &y, &z;

    __device__ inline __ref_float3(float& x, float& y, float& z)
	: x(x), y(y), z(z) {}

    __device__ inline __ref_float3& operator=(float3 const& v)
    {
	x = v.x;
	y = v.y;
	z = v.z;
	return *this;
    }

    __device__ inline __ref_float3& operator=(float4 const& v)
    {
	x = v.x;
	y = v.y;
	z = v.z;
	return *this;
    }

    __device__ inline operator float4() const
    {
	return make_float4(x, y, z, 0);
    }

    __device__ inline operator float3() const
    {
	return make_float3(x, y, z);
    }
};

struct __cref_float3
{
    float const &x, &y, &z;

    __device__ inline __cref_float3(float const& x, float const& y, float const& z)
	: x(x), y(y), z(z) {}

    __device__ inline operator float4() const
    {
	return make_float4(x, y, z, 0);
    }

    __device__ inline operator float3() const
    {
	return make_float3(x, y, z);
    }
};

struct __ref_float3_int
{
    float &x, &y, &z;
    int &w;

    __device__ inline __ref_float3_int(float& x, float& y, float& z, int& w)
	: x(x), y(y), z(z), w(w) {}

    __device__ inline __ref_float3_int& operator=(float4 const& v)
    {
	x = v.x;
	y = v.y;
	z = v.z;
	w = __float_as_int(v.w);
	return *this;
    }

    __device__ inline operator float4() const
    {
	return make_float4(x, y, z, __int_as_float(w));
    }
};

struct __cref_float3_int
{
    float const &x, &y, &z;
    int const &w;

    __device__ inline __cref_float3_int(float const& x, float const& y, float const& z, int const& w)
	: x(x), y(y), z(z), w(w) {}

    __device__ inline operator float4() const
    {
	return make_float4(x, y, z, __int_as_float(w));
    }
};

struct __ref_float4
{
    float &x, &y, &z, &w;

    __device__ inline __ref_float4(float& x, float& y, float& z, float& w)
	: x(x), y(y), z(z), w(w) {}

    __device__ inline __ref_float4& operator=(float4 const& v)
    {
	x = v.x;
	y = v.y;
	z = v.z;
	w = v.w;
	return *this;
    }

    __device__ inline operator float4() const
    {
	return make_float4(x, y, z, w);
    }
};

struct __cref_float4
{
    float const &x, &y, &z, &w;

    __device__ inline __cref_float4(float const& x, float const& y, float const& z, float const& w)
	: x(x), y(y), z(z), w(w) {}

    __device__ inline operator float4() const
    {
	return make_float4(x, y, z, w);
    }
};

__device__ inline __ref_float2_int operator,(vector<float, 2>& v, int& s)
{
    return __ref_float2_int(v.x, v.y, s);
}

__device__ inline __cref_float2_int operator,(vector<float, 2> const& v, int const& s)
{
    return __cref_float2_int(v.x, v.y, s);
}

__device__ inline __ref_float3 operator,(vector<float, 2>& v, float& s)
{
    return __ref_float3(v.x, v.y, s);
}

__device__ inline __cref_float3 operator,(vector<float, 2> const& v, float const& s)
{
    return __cref_float3(v.x, v.y, s);
}

__device__ inline __ref_float3_int operator,(vector<float, 3>& v, int& s)
{
    return __ref_float3_int(v.x, v.y, v.z, s);
}

__device__ inline __cref_float3_int operator,(vector<float, 3> const& v, int const& s)
{
    return __cref_float3_int(v.x, v.y, v.z, s);
}

__device__ inline __ref_float4 operator,(vector<float, 3>& v, float& s)
{
    return __ref_float4(v.x, v.y, v.z, s);
}

__device__ inline __cref_float4 operator,(vector<float, 3> const& v, float const& s)
{
    return __cref_float4(v.x, v.y, v.z, s);
}

}} // namespace ljgpu::cu

#endif /* ! LJGPU_MATH_GPU_PACK_CUH */
