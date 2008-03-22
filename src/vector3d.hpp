/* Three-dimensional vector
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

#ifndef _VECTOR3D_HPP
#define _VECTOR3D_HPP

#include <cmath>


/**
 * Three-dimensional vector
 */
template <typename T>
class vector3d
{
public:
    T x, y, z;

    /**
     * default constructor
     */
    vector3d()
    {
    }

    /**
     * copy constructor
     */
    vector3d(const vector3d<T>& v) : x(v.x), y(v.y), z(v.z)
    {
    }

    /**
     * value constructor
     */
    vector3d(const T& x, const T& y, const T& z) : x(x), y(y), z(z)
    {
    }

    size_t size() const
    {
	return 3;
    }

    /**
     * assignment operator
     */
    vector3d<T>& operator=(const vector3d<T>& v)
    {
	x = v.x;
	y = v.y;
	z = v.z;
	return *this;
    }

    vector3d<T>& operator=(const T& s)
    {
	x = s;
	y = s;
	z = s;
	return *this;
    }

    vector3d<T>& operator+=(const vector3d<T>& v)
    {
	x += v.x;
	y += v.y;
	z += v.z;
	return *this;
    }

    vector3d<T>& operator-=(const vector3d<T>& v)
    {
	x -= v.x;
	y -= v.y;
	z -= v.z;
	return *this;
    }

    vector3d<T>& operator*=(const T& s)
    {
	x *= s;
	y *= s;
	z *= s;
	return *this;
    }

    vector3d<T>& operator/=(const T& s)
    {
	x /= s;
	y /= s;
	z /= s;
	return *this;
    }

    vector3d<T> operator+(const vector3d<T>& v) const
    {
	return vector3d<T>(x + v.x, y + v.y, z + v.z);
    }

    vector3d<T> operator-(const vector3d<T>& v) const
    {
	return vector3d<T>(x - v.x, y - v.y, z - v.z);
    }

    vector3d<T> operator*(const T& s) const
    {
	return vector3d<T>(x * s, y * s, z * s);
    }

    vector3d<T> operator/(const T& s) const
    {
	return vector3d<T>(x / s, y / s, z / s);
    }

    /**
     * dot product
     */
    T operator*(const vector3d<T>& v) const
    {
	return x * v.x + y * v.y + z * v.z;
    }
};


/**
 * round vector components to nearest integer
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
 * round vector components down to nearest integer
 */
template <typename T>
vector3d<T> floor(const vector3d<T>& v)
{   
    return vector3d<T>(std::floor(v.x), std::floor(v.y), std::floor(v.z));
}


/**
 * round vector components up to nearest integer
 */
template <typename T>
vector3d<T> ceil(const vector3d<T>& v)
{   
    return vector3d<T>(std::ceil(v.x), std::ceil(v.y), std::ceil(v.z));
}


#endif /* ! _VECTOR3D_HPP */
