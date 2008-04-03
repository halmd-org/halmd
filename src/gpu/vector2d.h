/* 2-dimensional floating-point vector operations for CUDA device functions
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

#ifndef MDSIM_GPU_VECTOR2D_H
#define MDSIM_GPU_VECTOR2D_H

/**
 * componentwise vector addition
 */
__device__ float2 operator+(const float2& v, const float2& w)
{
    return make_float2(v.x + w.x, v.y + w.y);
}


/**
 * componentwise vector subtraction
 */
__device__ float2 operator-(const float2& v, const float2& w)
{
    return make_float2(v.x - w.x, v.y - w.y);
}


/**
 * scalar product
 */
__device__ float operator*(const float2& v, const float2& w)
{
    return v.x * w.x + v.y * w.y;
}


/**
 * scalar multiplication
 */
__device__ float2 operator*(const float2& v, const float& s)
{
    return make_float2(v.x * s, v.y * s);
}


/**
 * scalar multiplication
 */
__device__ float2 operator*(const float& s, const float2& v)
{
    return make_float2(s * v.x, s * v.y);
}


/**
 * scalar division
 */
__device__ float2 operator/(const float2& v, const float& s)
{
    return make_float2(v.x / s, v.y / s);
}


/**
 * assignment by componentwise vector addition
 */
__device__ float2& operator+=(float2& v, const float2& w)
{
    v.x += w.x;
    v.y += w.y;
    return v;
}


/**
 * assignment by componentwise vector subtraction
 */
__device__ float2& operator-=(float2& v, const float2& w)
{
    v.x -= w.x;
    v.y -= w.y;
    return v;
}


/**
 * assignment by scalar multiplication
 */
__device__ float2& operator*=(float2& v, const float& s)
{
    v.x *= s;
    v.y -= s;
    return v;
}


/**
 * assignment by scalar division
 */
__device__ float2& operator/=(float2& v, const float& s)
{
    v.x /= s;
    v.y /= s;
    return v;
}


/**
 * equality comparison operator
 */
__device__ bool operator==(const float2& v, const float2& w)
{
    return (v.x == w.x && v.y == w.y) ? true : false;
}


/**
 * returns vector with components set to given scalar
 */
template <typename T>
__device__ T make_floatn(const float& s);

template <>
__device__ float2 make_floatn(const float& s)
{
    return make_float2(s, s);
}


/**
 * componentwise round to nearest integer
 */
__device__ float2 rintf(const float2& v)
{
    return make_float2(rintf(v.x), rintf(v.y));
}


/**
 * componentwise round to nearest integer, away from zero
 */
__device__ float2 roundf(const float2& v)
{
    return make_float2(roundf(v.x), roundf(v.y));
}


/**
 * componentwise round to nearest integer not greater than argument
 */
__device__ float2 floorf(const float2& v)
{
    return make_float2(floorf(v.x), floorf(v.y));
}


/**
 * componentwise round to nearest integer not less argument
 */
__device__ float2 ceilf(const float2& v)
{
    return make_float2(ceilf(v.x), ceilf(v.y));
}


#endif /* ! MDSIM_GPU_VECTOR2D_H */
