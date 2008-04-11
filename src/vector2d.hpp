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
#include <iostream>


/**
 * 2-dimensional floating-point vector
 */
template <typename T>
class vector2d
{
public:
    T x, y;

public:
    vector2d()
    {
    }

    /**
     * initialization by vector
     */
    vector2d(const vector2d<T>& v) : x(v.x), y(v.y)
    {
    }

    /**
     * initialization by scalar
     */
    vector2d(const T& s) : x(s), y(s)
    {
    }

    /**
     * initialization by scalar components
     */
    vector2d(const T& x, const T& y) : x(x), y(y)
    {
    }

    /**
     * dimension of vector space
     */
    static unsigned int dim()
    {
	return 2;
    }

    /**
     * equality comparison
     */
    bool operator==(const vector2d<T>& v)
    {
	return (v.x == x && v.y == y);
    }

    /**
     * inequality comparison
     */
    bool operator!=(const vector2d<T>& v)
    {
	return (v.x != x || v.y != y);
    }

    /**
     * assignment by vector
     */
    vector2d<T>& operator=(const vector2d<T>& v)
    {
	x = v.x;
	y = v.y;
	return *this;
    }

    /**
     * assignment by scalar
     */
    vector2d<T>& operator=(const T& s)
    {
	x = s;
	y = s;
	return *this;
    }

    /**
     * assignment by componentwise vector addition
     */
    vector2d<T>& operator+=(const vector2d<T>& v)
    {
	x += v.x;
	y += v.y;
	return *this;
    }

    /**
     * assignment by componentwise vector subtraction
     */
    vector2d<T>& operator-=(const vector2d<T>& v)
    {
	x -= v.x;
	y -= v.y;
	return *this;
    }

    /**
     * assignment by scalar multiplication
     */
    vector2d<T>& operator*=(const T& s)
    {
	x *= s;
	y *= s;
	return *this;
    }

    /**
     * assignment by scalar division
     */
    vector2d<T>& operator/=(const T& s)
    {
	x /= s;
	y /= s;
	return *this;
    }

    /**
     * componentwise vector addition
     */
    vector2d<T> operator+(const vector2d<T>& v) const
    {
	return vector2d<T>(x + v.x, y + v.y);
    }

    /**
     * componentwise vector subtraction
     */
    vector2d<T> operator-(const vector2d<T>& v) const
    {
	return vector2d<T>(x - v.x, y - v.y);
    }

    /**
     * scalar product
     */
    T operator*(const vector2d<T>& v) const
    {
	return x * v.x + y * v.y;
    }

    /**
     * scalar multiplication
     */
    vector2d<T> operator*(const T& s) const
    {
	return vector2d<T>(x * s, y * s);
    }

    /**
     * scalar multiplication
     */
    friend vector2d<T> operator*(const T& s, const vector2d<T>& v)
    {
	return vector2d<T>(s * v.x, s * v.y);
    }

    /**
     * scalar division
     */
    vector2d<T> operator/(const T& s) const
    {
	return vector2d<T>(x / s, y / s);
    }

    /**
     * write vector components to output stream
     */
    friend std::ostream& operator<<(std::ostream& os, const vector2d<T>& v)
    {
	os << v.x << "\t" << v.y;
	return os;
    }

    /**
     * read vector components from input stream
     */
    friend std::istream& operator>>(std::istream& is, vector2d<T>& v)
    {
	is >> v.x >> v.y >> v.z;
	return is;
    }
};


/**
 * componentwise round to nearest integer
 */
template <typename T>
vector2d<T> rint(const vector2d<T>& v);

template <>
vector2d<float> rint(const vector2d<float>& v)
{
    return vector2d<float>(rintf(v.x), rintf(v.y));
}

template <>
vector2d<double> rint(const vector2d<double>& v)
{
    return vector2d<double>(rint(v.x), rint(v.y));
}


/**
 * componentwise round to nearest integer, away from zero
 */
template <typename T>
vector2d<T> round(const vector2d<T>& v);

template <>
vector2d<float> round(const vector2d<float>& v)
{
    return vector2d<float>(roundf(v.x), roundf(v.y));
}

template <>
vector2d<double> round(const vector2d<double>& v)
{
    return vector2d<double>(round(v.x), round(v.y));
}


/**
 * componentwise round to nearest integer not greater than argument
 */
template <typename T>
vector2d<T> floor(const vector2d<T>& v);

template <>
vector2d<float> floor(const vector2d<float>& v)
{
    return vector2d<float>(floorf(v.x), floorf(v.y));
}

template <>
vector2d<double> floor(const vector2d<double>& v)
{
    return vector2d<double>(floor(v.x), floor(v.y));
}


/**
 * componentwise round to nearest integer not less argument
 */
template <typename T>
vector2d<T> ceil(const vector2d<T>& v);

template <>
vector2d<float> ceil(const vector2d<float>& v)
{
    return vector2d<float>(ceilf(v.x), ceilf(v.y));
}

template <>
vector2d<double> ceil(const vector2d<double>& v)
{
    return vector2d<double>(ceil(v.x), ceil(v.y));
}


/**
 * componentwise square root function
 */
template <typename T>
vector2d<T> sqrt(const vector2d<T>& v);

template <>
vector2d<float> sqrt(const vector2d<float>& v)
{
    return vector2d<float>(sqrtf(v.x), sqrtf(v.y));
}

template <>
vector2d<double> sqrt(const vector2d<double>& v)
{
    return vector2d<double>(sqrt(v.x), sqrt(v.y));
}

#endif /* ! MDSIM_VECTOR2D_HPP */
