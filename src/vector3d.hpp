/* 3-dimensional floating-point vector
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

#ifndef MDSIM_VECTOR3D_HPP
#define MDSIM_VECTOR3D_HPP

#include <math.h>
#include <cuda/vector_types.h>
#include <cuda/vector_functions.h>


/**
 * 3-dimensional floating-point vector
 */
class vector3d
{
public:
    float x, y, z;

public:
    vector3d()
    {
    }

    /**
     * initialization by vector
     */
    vector3d(const vector3d& v) : x(v.x), y(v.y), z(v.z)
    {
    }

    /**
     * initialization by vector
     */
    vector3d(const float3& v) : x(v.x), y(v.y), z(v.z)
    {
    }

    /**
     * initialization by scalar
     */
    vector3d(const float& s) : x(s), y(s), z(s)
    {
    }

    /**
     * initialization by scalar components
     */
    vector3d(const float& x, const float& y, const float& z) : x(x), y(y), z(z)
    {
    }

    /**
     * conversion to CUDA vector type
     */
    float3 floatn() const
    {
	return make_float3(x, y, z);
    }

    /**
     * dimension of vector space
     */
    unsigned int dim() const
    {
	return 3;
    }

    /**
     * assignment by vector
     */
    vector3d& operator=(const vector3d& v)
    {
	x = v.x;
	y = v.y;
	z = v.z;
	return *this;
    }

    /**
     * assignment by vector
     */
    vector3d& operator=(const float3& v)
    {
	x = v.x;
	y = v.y;
	z = v.z;
	return *this;
    }

    /**
     * assignment by scalar
     */
    vector3d& operator=(const float& s)
    {
	x = s;
	y = s;
	z = s;
	return *this;
    }

    /**
     * assignment by componentwise vector addition
     */
    vector3d& operator+=(const vector3d& v)
    {
	x += v.x;
	y += v.y;
	z += v.z;
	return *this;
    }

    /**
     * assignment by componentwise vector addition
     */
    vector3d& operator+=(const float3& v)
    {
	x += v.x;
	y += v.y;
	z += v.z;
	return *this;
    }

    /**
     * assignment by componentwise vector subtraction
     */
    vector3d& operator-=(const vector3d& v)
    {
	x -= v.x;
	y -= v.y;
	z -= v.z;
	return *this;
    }

    /**
     * assignment by componentwise vector subtraction
     */
    vector3d& operator-=(const float3& v)
    {
	x -= v.x;
	y -= v.y;
	z -= v.z;
	return *this;
    }

    /**
     * assignment by scalar multiplication
     */
    vector3d& operator*=(const float& s)
    {
	x *= s;
	y *= s;
	z *= s;
	return *this;
    }

    /**
     * assignment by scalar division
     */
    vector3d& operator/=(const float& s)
    {
	x /= s;
	y /= s;
	z /= s;
	return *this;
    }

    /**
     * componentwise vector addition
     */
    vector3d operator+(const vector3d& v) const
    {
	return vector3d(x + v.x, y + v.y, z + v.z);
    }

    /**
     * componentwise vector subtraction
     */
    vector3d operator-(const vector3d& v) const
    {
	return vector3d(x - v.x, y - v.y, z - v.z);
    }

    /**
     * scalar product
     */
    float operator*(const vector3d& v) const
    {
	return x * v.x + y * v.y + z * v.z;
    }

    /**
     * scalar multiplication
     */
    vector3d operator*(const float& s) const
    {
	return vector3d(x * s, y * s, z * s);
    }

    /**
     * scalar division
     */
    vector3d operator/(const float& s) const
    {
	return vector3d(x / s, y / s, z / s);
    }

    /**
     * componentwise round to nearest integer
     */
    friend vector3d rintf(const vector3d& v)
    {
	return vector3d(rintf(v.x), rintf(v.y), rintf(v.z));
    }

    /*
     * componentwise round to nearest integer, away from zero
     */
    friend vector3d roundf(const vector3d& v)
    {
	return vector3d(roundf(v.x), roundf(v.y), roundf(v.z));
    }

    /**
     * componentwise round to nearest integer not greater than argument
     */
    friend vector3d floorf(const vector3d& v)
    {   
	return vector3d(floorf(v.x), floorf(v.y), floorf(v.z));
    }

    /**
     * componentwise round to nearest integer not less argument
     */
    friend vector3d ceilf(const vector3d& v)
    {   
	return vector3d(ceilf(v.x), ceilf(v.y), ceilf(v.z));
    }
};

#endif /* ! MDSIM_VECTOR3D_HPP */
