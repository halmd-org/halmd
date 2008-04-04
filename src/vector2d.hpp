/* 2-dimensional floating-point vector
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

#ifndef MDSIM_VECTOR2D_HPP
#define MDSIM_VECTOR2D_HPP

#include <math.h>
#include <cuda/vector_types.h>
#include <cuda/vector_functions.h>


/**
 * 2-dimensional floating-point vector
 */
class vector2d
{
public:
    float x, y;

public:
    vector2d()
    {
    }

    /**
     * initialization by vector
     */
    vector2d(const vector2d& v) : x(v.x), y(v.y)
    {
    }

    /**
     * initialization by vector
     */
    vector2d(const float2& v) : x(v.x), y(v.y)
    {
    }

    /**
     * initialization by scalar
     */
    vector2d(const float& s) : x(s), y(s)
    {
    }

    /**
     * initialization by scalar components
     */
    vector2d(const float& x, const float& y) : x(x), y(y)
    {
    }

    /**
     * conversion to CUDA vector type
     */
    float2 floatn() const
    {
	return make_float2(x, y);
    }

    /**
     * dimension of vector space
     */
    unsigned int dim() const
    {
	return 2;
    }

    /**
     * assignment by vector
     */
    vector2d& operator=(const vector2d& v)
    {
	x = v.x;
	y = v.y;
	return *this;
    }

    /**
     * assignment by vector
     */
    vector2d& operator=(const float2& v)
    {
	x = v.x;
	y = v.y;
	return *this;
    }

    /**
     * assignment by scalar
     */
    vector2d& operator=(const float& s)
    {
	x = s;
	y = s;
	return *this;
    }

    /**
     * assignment by componentwise vector addition
     */
    vector2d& operator+=(const vector2d& v)
    {
	x += v.x;
	y += v.y;
	return *this;
    }

    /**
     * assignment by componentwise vector addition
     */
    vector2d& operator+=(float2& v)
    {
	x += v.x;
	y += v.y;
	return *this;
    }

    /**
     * assignment by componentwise vector subtraction
     */
    vector2d& operator-=(const vector2d& v)
    {
	x -= v.x;
	y -= v.y;
	return *this;
    }

    /**
     * assignment by componentwise vector subtraction
     */
    vector2d& operator-=(float2& v)
    {
	x -= v.x;
	y -= v.y;
	return *this;
    }

    /**
     * assignment by scalar multiplication
     */
    vector2d& operator*=(const float& s)
    {
	x *= s;
	y *= s;
	return *this;
    }

    /**
     * assignment by scalar division
     */
    vector2d& operator/=(const float& s)
    {
	x /= s;
	y /= s;
	return *this;
    }

    /**
     * componentwise vector addition
     */
    vector2d operator+(const vector2d& v) const
    {
	return vector2d(x + v.x, y + v.y);
    }

    /**
     * componentwise vector subtraction
     */
    vector2d operator-(const vector2d& v) const
    {
	return vector2d(x - v.x, y - v.y);
    }

    /**
     * scalar product
     */
    float operator*(const vector2d& v) const
    {
	return x * v.x + y * v.y;
    }

    /**
     * scalar multiplication
     */
    vector2d operator*(const float& s) const
    {
	return vector2d(x * s, y * s);
    }

    /**
     * scalar division
     */
    vector2d operator/(const float& s) const
    {
	return vector2d(x / s, y / s);
    }

    /**
     * componentwise round to nearest integer
     */
    friend vector2d rintf(const vector2d& v)
    {
	return vector2d(rintf(v.x), rintf(v.y));
    }

    /**
     * componentwise round to nearest integer, away from zero
     */
    friend vector2d roundf(const vector2d& v)
    {
	return vector2d(roundf(v.x), roundf(v.y));
    }

    /**
     * componentwise round to nearest integer not greater than argument
     */
    friend vector2d floorf(const vector2d& v)
    {
	return vector2d(floorf(v.x), floorf(v.y));
    }

    /**
     * componentwise round to nearest integer not less argument
     */
    friend vector2d ceilf(const vector2d& v)
    {
	return vector2d(ceilf(v.x), ceilf(v.y));
    }
};

#endif /* ! MDSIM_VECTOR2D_HPP */
