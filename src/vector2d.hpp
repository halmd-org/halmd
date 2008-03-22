/* Two-dimensional vector
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

#ifndef _VECTOR2D_HPP
#define _VECTOR2D_HPP

#include <cmath>


/**
 * Two-dimensional vector
 */
template <typename T>
class vector2d
{
public:
    T x, y;

    /**
     * default constructor
     */
    vector2d()
    {
    }

    /**
     * copy constructor
     */
    vector2d(const vector2d<T>& v) : x(v.x), y(v.y)
    {
    }

    /**
     * value constructor
     */
    vector2d(const T& x, const T& y) : x(x), y(y)
    {
    }

    size_t size() const
    {
	return 2;
    }

    /**
     * assignment operator
     */
    vector2d<T>& operator=(const vector2d<T>& v)
    {
	x = v.x;
	y = v.y;
	return *this;
    }

    vector2d<T>& operator=(const T& s)
    {
	x = s;
	y = s;
	return *this;
    }

    vector2d<T>& operator+=(const vector2d<T>& v)
    {
	x += v.x;
	y += v.y;
	return *this;
    }

    vector2d<T>& operator-=(const vector2d<T>& v)
    {
	x -= v.x;
	y -= v.y;
	return *this;
    }

    vector2d<T>& operator*=(const T& s)
    {
	x *= s;
	y *= s;
	return *this;
    }

    vector2d<T>& operator/=(const T& s)
    {
	x /= s;
	y /= s;
	return *this;
    }

    vector2d<T> operator+(const vector2d<T>& v) const
    {
	return vector2d<T>(x + v.x, y + v.y);
    }

    vector2d<T> operator-(const vector2d<T>& v) const
    {
	return vector2d<T>(x - v.x, y - v.y);
    }

    vector2d<T> operator*(const T& s) const
    {
	return vector2d<T>(x * s, y * s);
    }

    vector2d<T> operator/(const T& s) const
    {
	return vector2d<T>(x / s, y / s);
    }

    /**
     * dot product
     */
    T operator*(const vector2d<T>& v) const
    {
	return x * v.x + y * v.y;
    }
};


/**
 * round vector components to nearest integer
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
 * round vector components down to nearest integer
 */
template <typename T>
vector2d<T> floor(const vector2d<T>& v)
{
    return vector2d<T>(std::floor(v.x), std::floor(v.y));
}


/**
 * round vector components up to nearest integer
 */
template <typename T>
vector2d<T> ceil(const vector2d<T>& v)
{
    return vector2d<T>(std::ceil(v.x), std::ceil(v.y));
}


#endif /* ! _VECTOR2D_HPP */
