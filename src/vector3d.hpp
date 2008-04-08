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
#include <iostream>


/**
 * 3-dimensional floating-point vector
 */
template <typename T>
class vector3d
{
public:
    T x, y, z;

public:
    vector3d()
    {
    }

    /**
     * initialization by vector
     */
    vector3d(const vector3d<T>& v) : x(v.x), y(v.y), z(v.z)
    {
    }

    /**
     * initialization by scalar
     */
    vector3d(const T& s) : x(s), y(s), z(s)
    {
    }

    /**
     * initialization by scalar components
     */
    vector3d(const T& x, const T& y, const T& z) : x(x), y(y), z(z)
    {
    }

    /**
     * dimension of vector space
     */
    static unsigned int dim()
    {
	return 3;
    }

    /**
     * assignment by vector
     */
    vector3d<T>& operator=(const vector3d<T>& v)
    {
	x = v.x;
	y = v.y;
	z = v.z;
	return *this;
    }

    /**
     * assignment by scalar
     */
    vector3d<T>& operator=(const T& s)
    {
	x = s;
	y = s;
	z = s;
	return *this;
    }

    /**
     * assignment by componentwise vector addition
     */
    vector3d<T>& operator+=(const vector3d<T>& v)
    {
	x += v.x;
	y += v.y;
	z += v.z;
	return *this;
    }

    /**
     * assignment by componentwise vector subtraction
     */
    vector3d<T>& operator-=(const vector3d<T>& v)
    {
	x -= v.x;
	y -= v.y;
	z -= v.z;
	return *this;
    }

    /**
     * assignment by scalar multiplication
     */
    vector3d<T>& operator*=(const T& s)
    {
	x *= s;
	y *= s;
	z *= s;
	return *this;
    }

    /**
     * assignment by scalar division
     */
    vector3d<T>& operator/=(const T& s)
    {
	x /= s;
	y /= s;
	z /= s;
	return *this;
    }

    /**
     * componentwise vector addition
     */
    vector3d<T> operator+(const vector3d<T>& v) const
    {
	return vector3d<T>(x + v.x, y + v.y, z + v.z);
    }

    /**
     * componentwise vector subtraction
     */
    vector3d<T> operator-(const vector3d<T>& v) const
    {
	return vector3d<T>(x - v.x, y - v.y, z - v.z);
    }

    /**
     * scalar product
     */
    T operator*(const vector3d<T>& v) const
    {
	return x * v.x + y * v.y + z * v.z;
    }

    /**
     * scalar multiplication
     */
    vector3d<T> operator*(const T& s) const
    {
	return vector3d<T>(x * s, y * s, z * s);
    }

    /**
     * scalar multiplication
     */
    friend vector3d<T> operator*(const T& s, const vector3d<T>& v)
    {
	return vector3d<T>(s * v.x, s * v.y, s * v.z);
    }

    /**
     * scalar division
     */
    vector3d<T> operator/(const T& s) const
    {
	return vector3d<T>(x / s, y / s, z / s);
    }

    /**
     * write vector components to output stream
     */
    friend std::ostream& operator<<(std::ostream& os, const vector3d<T>& v)
    {
	os << v.x << "\t" << v.y << "\t" << v.z;
	return os;
    }

    /**
     * read vector components from input stream
     */
    friend std::istream& operator>>(std::istream& is, vector3d<T>& v)
    {
	is >> v.x >> v.y >> v.z;
	return is;
    }
};


/**
 * componentwise round to nearest integer
 */
template <typename T>
vector3d<T> rint(const vector3d<T>& v);

template <>
vector3d<float> rint(const vector3d<float>& v)
{
    return vector3d<float>(rintf(v.x), rintf(v.y), rintf(v.z));
}

template <>
vector3d<double> rint(const vector3d<double>& v)
{
    return vector3d<double>(rint(v.x), rint(v.y), rint(v.z));
}


/**
 * componentwise round to nearest integer, away from zero
 */
template <typename T>
vector3d<T> round(const vector3d<T>& v);

template <>
vector3d<float> round(const vector3d<float>& v)
{
    return vector3d<float>(roundf(v.x), roundf(v.y), roundf(v.z));
}

template <>
vector3d<double> round(const vector3d<double>& v)
{
    return vector3d<double>(round(v.x), round(v.y), round(v.z));
}


/**
 * componentwise round to nearest integer not greater than argument
 */
template <typename T>
vector3d<T> floor(const vector3d<T>& v);

template <>
vector3d<float> floor(const vector3d<float>& v)
{   
    return vector3d<float>(floorf(v.x), floorf(v.y), floorf(v.z));
}

template <>
vector3d<double> floor(const vector3d<double>& v)
{   
    return vector3d<double>(floor(v.x), floor(v.y), floor(v.z));
}


/**
 * componentwise round to nearest integer not less argument
 */
template <typename T>
vector3d<T> ceil(const vector3d<T>& v);

template <>
vector3d<float> ceil(const vector3d<float>& v)
{   
    return vector3d<float>(ceilf(v.x), ceilf(v.y), ceilf(v.z));
}

template <>
vector3d<double> ceil(const vector3d<double>& v)
{   
    return vector3d<double>(ceil(v.x), ceil(v.y), ceil(v.z));
}


/**
 * componentwise square root function
 */
template <typename T>
vector3d<T> sqrt(const vector3d<T>& v);

template <>
vector3d<float> sqrt(const vector3d<float>& v)
{   
    return vector3d<float>(sqrtf(v.x), sqrtf(v.y), sqrtf(v.z));
}

template <>
vector3d<double> sqrt(const vector3d<double>& v)
{   
    return vector3d<double>(sqrt(v.x), sqrt(v.y), sqrt(v.z));
}

#endif /* ! MDSIM_VECTOR3D_HPP */
