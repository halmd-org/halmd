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

template <typename T1, typename T2>
class __pack;

template <>
class __pack<vector<float, 2>&, int&>
{
public:
    typedef vector<float, 2> vector_type;
    typedef __pack<vector_type&, int&> type;

    __device__ inline __pack(vector_type& v, int& s): v(v), s(s) {}

    __device__ inline type& operator=(float3 const& w)
    {
	v = vector_type(w.x, w.y);
	s = __float_as_int(w.z);
	return *this;
    }

    __device__ inline type& operator=(float4 const& w)
    {
	v = vector_type(w.x, w.y);
	s = __float_as_int(w.z);
	return *this;
    }

    __device__ inline operator float4() const
    {
	return make_float4(v.x, v.y, __int_as_float(s), 0);
    }

    __device__ inline operator float3() const
    {
	return make_float3(v.x, v.y, __int_as_float(s));
    }

private:
    vector_type& v;
    int& s;
};

template <>
class __pack<vector<float, 2> const&, int const&>
{
public:
    typedef vector<float, 2> vector_type;
    typedef __pack<vector_type const&, int const&> type;

    __device__ inline __pack(vector_type const& v, int const& s): v(v), s(s) {}

    __device__ inline operator float4() const
    {
	return make_float4(v.x, v.y, __int_as_float(s), 0);
    }

    __device__ inline operator float3() const
    {
	return make_float3(v.x, v.y, __int_as_float(s));
    }

private:
    vector_type const& v;
    int const& s;
};

template <>
class __pack<vector<float, 3>&, int&>
{
public:
    typedef vector<float, 3> vector_type;
    typedef __pack<vector_type&, int&> type;

    __device__ inline __pack(vector_type& v, int& s): v(v), s(s) {}

    __device__ inline type& operator=(float4 const& w)
    {
	v = vector_type(w.x, w.y, w.z);
	s = __float_as_int(w.w);
	return *this;
    }

    __device__ inline operator float4() const
    {
	return make_float4(v.x, v.y, v.z, __int_as_float(s));
    }

private:
    vector_type& v;
    int& s;
};

template <>
class __pack<vector<float, 3> const&, int const&>
{
public:
    typedef vector<float, 3> vector_type;
    typedef __pack<vector_type const&, int const&> type;

    __device__ inline __pack(vector_type const& v, int const& s): v(v), s(s) {}

    __device__ inline operator float4() const
    {
	return make_float4(v.x, v.y, v.z, __int_as_float(s));
    }

private:
    vector_type const& v;
    int const& s;
};

__device__ inline __pack<vector<float, 2>&, int&>
operator,(vector<float, 2>& v, int& s)
{
    return __pack<vector<float, 2>&, int&>(v, s);
}

__device__ inline __pack<vector<float, 2> const&, int const&>
operator,(vector<float, 2> const& v, int const& s)
{
    return __pack<vector<float, 2> const&, int const&>(v, s);
}

__device__ inline __pack<vector<float, 3>&, int&>
operator,(vector<float, 3>& v, int& s)
{
    return __pack<vector<float, 3>&, int&>(v, s);
}

__device__ inline __pack<vector<float, 3> const&, int const&>
operator,(vector<float, 3> const& v, int const& s)
{
    return __pack<vector<float, 3> const&, int const&>(v, s);
}

}} // namespace ljgpu::cu

#endif /* ! LJGPU_MATH_GPU_PACK_CUH */
