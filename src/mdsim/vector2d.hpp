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

#include <boost/array.hpp>
#include <cmath>
#include <iostream>
#include <cuda_wrapper.hpp>
#include <xdr/iostream.hpp>


/**
 * 2-dimensional floating-point vector
 */
template <typename T>
class vector2d : public boost::array<T, 2>
{
public:
    typedef T value_type;

public:
    vector2d()
    {
    }

    /**
     * initialization by vector
     */
    vector2d(vector2d<T> const& v)
    {
	(*this)[0] = v[0];
	(*this)[1] = v[1];
    }

    /**
     * initialization by coalesced GPU floating-point vector
     */
    vector2d(float2 const& v)
    {
	(*this)[0] = v.x;
	(*this)[1] = v.y;
    }

    /**
     * initialization by scalar
     */
    vector2d(T const& s)
    {
	(*this)[0] = s;
	(*this)[1] = s;
    }

    /**
     * initialization by scalar components
     */
    vector2d(T const& x, T const& y)
    {
	(*this)[0] = x;
	(*this)[1] = y;
    }

    /**
     * equality comparison
     */
    bool operator==(vector2d<T> const& v) const
    {
	return (v[0] == (*this)[0] && v[1] == (*this)[1]);
    }

    /**
     * inequality comparison
     */
    bool operator!=(vector2d<T> const& v) const
    {
	return (v[0] != (*this)[0] || v[1] != (*this)[1]);
    }

    /**
     * componentwise less than comparison
     */
    bool operator<(vector2d<T> const& v) const
    {
	return (v[0] < (*this)[0] && v[1] < (*this)[1]);
    }

    /**
     * componentwise greater than comparison
     */
    bool operator>(vector2d<T> const& v) const
    {
	return (v[0] > (*this)[0] && v[1] > (*this)[1]);
    }

    /**
     * componentwise less than or equal to comparison
     */
    bool operator<=(vector2d<T> const& v) const
    {
	return (v[0] <= (*this)[0] && v[1] <= (*this)[1]);
    }

    /**
     * componentwise greater than or equal to comparison
     */
    bool operator>=(vector2d<T> const& v) const
    {
	return (v[0] >= (*this)[0] && v[1] >= (*this)[1]);
    }

    /**
     * assignment by vector
     */
    vector2d<T>& operator=(vector2d<T> const& v)
    {
	(*this)[0] = v[0];
	(*this)[1] = v[1];
	return *this;
    }

    /**
     * assignment by scalar
     */
    vector2d<T>& operator=(T const& s)
    {
	(*this)[0] = s;
	(*this)[1] = s;
	return *this;
    }

    /**
     * assignment by componentwise vector addition
     */
    vector2d<T>& operator+=(vector2d<T> const& v)
    {
	(*this)[0] += v[0];
	(*this)[1] += v[1];
	return *this;
    }

    /**
     * assignment by componentwise vector subtraction
     */
    vector2d<T>& operator-=(vector2d<T> const& v)
    {
	(*this)[0] -= v[0];
	(*this)[1] -= v[1];
	return *this;
    }

    /**
     * assignment by scalar multiplication
     */
    vector2d<T>& operator*=(T const& s)
    {
	(*this)[0] *= s;
	(*this)[1] *= s;
	return *this;
    }

    /**
     * assignment by scalar division
     */
    vector2d<T>& operator/=(T const& s)
    {
	(*this)[0] /= s;
	(*this)[1] /= s;
	return *this;
    }

    /**
     * componentwise vector addition
     */
    friend vector2d<T> operator+(vector2d<T> v, vector2d<T> const& w)
    {
	v[0] += w[0];
	v[1] += w[1];
	return v;
    }

    /**
     * componentwise vector subtraction
     */
    friend vector2d<T> operator-(vector2d<T> v, vector2d<T> const& w)
    {
	v[0] -= w[0];
	v[1] -= w[1];
	return v;
    }

    /**
     * scalar product
     */
    T operator*(vector2d<T> const& v) const
    {
	return (*this)[0] * v[0] + (*this)[1] * v[1];
    }

    /**
     * scalar multiplication
     */
    friend vector2d<T> operator*(vector2d<T> v, T const& s)
    {
	v[0] *= s;
	v[1] *= s;
	return v;
    }

    /**
     * scalar multiplication
     */
    friend vector2d<T> operator*(T const& s, vector2d<T> v)
    {
	v[0] *= s;
	v[1] *= s;
	return v;
    }

    /**
     * scalar division
     */
    friend vector2d<T> operator/(vector2d<T> v, T const& s)
    {
	v[0] /= s;
	v[1] /= s;
	return v;
    }

    /**
     * write vector components to output stream
     */
    friend std::ostream& operator<<(std::ostream& os, vector2d<T> const& v)
    {
	os << v[0] << "\t" << v[1];
	return os;
    }

    /**
     * read vector components from input stream
     */
    friend std::istream& operator>>(std::istream& is, vector2d<T>& v)
    {
	is >> v[0] >> v[1];
	return is;
    }

    /**
     * translate from vector to XDR vector
     */
    friend xdr::ostream& operator<<(xdr::ostream& xdrs, vector2d<T> const& v)
    {
	xdrs << v[0] << v[1];
	return xdrs;
    }

    /**
     * translate to vector from XDR vector
     */
    friend xdr::istream& operator>>(xdr::istream& xdrs, vector2d<T>& v)
    {
	xdrs >> v[0] >> v[1];
	return xdrs;
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
    v[0] = rintf(v[0]);
    v[1] = rintf(v[1]);
    return v;
}

template <>
vector2d<double> rint(vector2d<double> v)
{
    v[0] = rint(v[0]);
    v[1] = rint(v[1]);
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
    v[0] = roundf(v[0]);
    v[1] = roundf(v[1]);
    return v;
}

template <>
vector2d<double> round(vector2d<double> v)
{
    v[0] = round(v[0]);
    v[1] = round(v[1]);
    return v;
}

/**
 * componentwise round to nearest integer not greater than argument
 */
template <typename T>
vector2d<T> floor(vector2d<T> v)
{
    v[0] = std::floor(v[0]);
    v[1] = std::floor(v[1]);
    return v;
}

/**
 * componentwise round to nearest integer not less argument
 */
template <typename T>
vector2d<T> ceil(vector2d<T> v)
{
    v[0] = std::ceil(v[0]);
    v[1] = std::ceil(v[1]);
    return v;
}

/**
 * componentwise square root function
 */
template <typename T>
vector2d<T> sqrt(vector2d<T> v)
{
    v[0] = std::sqrt(v[0]);
    v[1] = std::sqrt(v[1]);
    return v;
}

/**
 * componentwise cos function
 */
template <typename T>
vector2d<T> cos(vector2d<T> v)
{
    v[0] = std::cos(v[0]);
    v[1] = std::cos(v[1]);
    return v;
}

/**
 * componentwise sin function
 */
template <typename T>
vector2d<T> sin(vector2d<T> v)
{
    v[0] = std::sin(v[0]);
    v[1] = std::sin(v[1]);
    return v;
}

/**
 * convert to GPU floating-point type
 */
template <typename T>
float2 make_float(vector2d<T> const& v)
{
    return make_float2(v[0], v[1]);
}

#endif /* ! MDSIM_VECTOR2D_HPP */
