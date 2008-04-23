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
    vector2d(vector2d<T> const& v) : x(v.x), y(v.y)
    {
    }

    /**
     * initialization by scalar
     */
    vector2d(T const& s) : x(s), y(s)
    {
    }

    /**
     * initialization by scalar components
     */
    vector2d(T const& x, T const& y) : x(x), y(y)
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
    bool operator==(vector2d<T> const& v) const
    {
	return (v.x == x && v.y == y);
    }

    /**
     * inequality comparison
     */
    bool operator!=(vector2d<T> const& v) const
    {
	return (v.x != x || v.y != y);
    }

    /**
     * componentwise less than comparison
     */
    bool operator<(vector2d<T> const& v) const
    {
	return (v.x < x && v.y < y);
    }

    /**
     * componentwise greater than comparison
     */
    bool operator>(vector2d<T> const& v) const
    {
	return (v.x > x && v.y > y);
    }

    /**
     * componentwise less than or equal to comparison
     */
    bool operator<=(vector2d<T> const& v) const
    {
	return (v.x <= x && v.y <= y);
    }

    /**
     * componentwise greater than or equal to comparison
     */
    bool operator>=(vector2d<T> const& v) const
    {
	return (v.x >= x && v.y >= y);
    }

    /**
     * assignment by vector
     */
    vector2d<T>& operator=(vector2d<T> const& v)
    {
	x = v.x;
	y = v.y;
	return *this;
    }

    /**
     * assignment by scalar
     */
    vector2d<T>& operator=(T const& s)
    {
	x = s;
	y = s;
	return *this;
    }

    /**
     * assignment by componentwise vector addition
     */
    vector2d<T>& operator+=(vector2d<T> const& v)
    {
	x += v.x;
	y += v.y;
	return *this;
    }

    /**
     * assignment by componentwise vector subtraction
     */
    vector2d<T>& operator-=(vector2d<T> const& v)
    {
	x -= v.x;
	y -= v.y;
	return *this;
    }

    /**
     * assignment by scalar multiplication
     */
    vector2d<T>& operator*=(T const& s)
    {
	x *= s;
	y *= s;
	return *this;
    }

    /**
     * assignment by scalar division
     */
    vector2d<T>& operator/=(T const& s)
    {
	x /= s;
	y /= s;
	return *this;
    }

    /**
     * componentwise vector addition
     */
    friend vector2d<T> operator+(vector2d<T> v, vector2d<T> const& w)
    {
	v.x += w.x;
	v.y += w.y;
	return v;
    }

    /**
     * componentwise vector subtraction
     */
    friend vector2d<T> operator-(vector2d<T> v, vector2d<T> const& w)
    {
	v.x -= w.x;
	v.y -= w.y;
	return v;
    }

    /**
     * scalar product
     */
    T operator*(vector2d<T> const& v) const
    {
	return x * v.x + y * v.y;
    }

    /**
     * scalar multiplication
     */
    friend vector2d<T> operator*(vector2d<T> v, T const& s)
    {
	v.x *= s;
	v.y *= s;
	return v;
    }

    /**
     * scalar multiplication
     */
    friend vector2d<T> operator*(T const& s, vector2d<T> v)
    {
	v.x *= s;
	v.y *= s;
	return v;
    }

    /**
     * scalar division
     */
    friend vector2d<T> operator/(vector2d<T> v, T const& s)
    {
	v.x /= s;
	v.y /= s;
	return v;
    }

    /**
     * write vector components to output stream
     */
    friend std::ostream& operator<<(std::ostream& os, vector2d<T> const& v)
    {
	os << v.x << "\t" << v.y;
	return os;
    }

    /**
     * read vector components from input stream
     */
    friend std::istream& operator>>(std::istream& is, vector2d<T>& v)
    {
	is >> v.x >> v.y;
	return is;
    }
};


/**
 * componentwise round to nearest integer
 */
template <typename T>
vector2d<T> rint(vector2d<T> v);

template <>
vector2d<float> rint(vector2d<float> v)
{
    v.x = rintf(v.x);
    v.y = rintf(v.y);
    return v;
}

template <>
vector2d<double> rint(vector2d<double> v)
{
    v.x = rint(v.x);
    v.y = rint(v.y);
    return v;
}


/**
 * componentwise round to nearest integer, away from zero
 */
template <typename T>
vector2d<T> round(vector2d<T> v);

template <>
vector2d<float> round(vector2d<float> v)
{
    v.x = roundf(v.x);
    v.y = roundf(v.y);
    return v;
}

template <>
vector2d<double> round(vector2d<double> v)
{
    v.x = round(v.x);
    v.y = round(v.y);
    return v;
}


/**
 * componentwise round to nearest integer not greater than argument
 */
template <typename T>
vector2d<T> floor(vector2d<T> v);

template <>
vector2d<float> floor(vector2d<float> v)
{
    v.x = floorf(v.x);
    v.y = floorf(v.y);
    return v;
}

template <>
vector2d<double> floor(vector2d<double> v)
{
    v.x = floor(v.x);
    v.y = floor(v.y);
    return v;
}


/**
 * componentwise round to nearest integer not less argument
 */
template <typename T>
vector2d<T> ceil(vector2d<T> v);

template <>
vector2d<float> ceil(vector2d<float> v)
{
    v.x = ceilf(v.x);
    v.y = ceilf(v.y);
    return v;
}

template <>
vector2d<double> ceil(vector2d<double> v)
{
    v.x = ceil(v.x);
    v.y = ceil(v.y);
    return v;
}


/**
 * componentwise square root function
 */
template <typename T>
vector2d<T> sqrt(vector2d<T> v);

template <>
vector2d<float> sqrt(vector2d<float> v)
{
    v.x = sqrtf(v.x);
    v.y = sqrtf(v.y);
    return v;
}

template <>
vector2d<double> sqrt(vector2d<double> v)
{
    v.x = sqrt(v.x);
    v.y = sqrt(v.y);
    return v;
}

#endif /* ! MDSIM_VECTOR2D_HPP */
