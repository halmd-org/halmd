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

#include <boost/array.hpp>
#include <cmath>
#include <iostream>
#include <cuda_wrapper.hpp>
#include <xdr/iostream.hpp>


/**
 * 3-dimensional floating-point vector
 */
template <typename T>
class vector3d : public boost::array<T, 3>
{
public:
    typedef T value_type;

public:
    vector3d()
    {
    }

    /**
     * initialization by vector
     */
    vector3d(vector3d<T> const& v)
    {
	(*this)[0] = v[0];
	(*this)[1] = v[1];
	(*this)[2] = v[2];
    }

    /**
     * initialization by vector
     */
    vector3d(float4 const& v)
    {
	(*this)[0] = v.x;
	(*this)[1] = v.y;
	(*this)[2] = v.z;
    }

    /**
     * initialization by scalar
     */
    vector3d(T const& s)
    {
	(*this)[0] = s;
	(*this)[1] = s;
	(*this)[2] = s;
    }

    /**
     * initialization by scalar components
     */
    vector3d(T const& x, T const& y, T const& z)
    {
	(*this)[0] = x;
	(*this)[1] = y;
	(*this)[2] = z;
    }

    /**
     * equality comparison
     */
    bool operator==(vector3d<T> const& v) const
    {
        return (v[0] == (*this)[0] && v[1] == (*this)[1] && v[2] == (*this)[2]);
    }

    /**
     * inequality comparison
     */
    bool operator!=(vector3d<T> const& v) const
    {
        return (v[0] != (*this)[0] || v[1] != (*this)[1] || v[2] != (*this)[2]);
    }

    /**
     * componentwise less than comparison
     */
    bool operator<(vector3d<T> const& v) const
    {
	return (v[0] < (*this)[0] && v[1] < (*this)[1] && v[2] < (*this)[2]);
    }

    /**
     * componentwise greater than comparison
     */
    bool operator>(vector3d<T> const& v) const
    {
	return (v[0] > (*this)[0] && v[1] > (*this)[1] && v[2] > (*this)[2]);
    }

    /**
     * componentwise less than or equal to comparison
     */
    bool operator<=(vector3d<T> const& v) const
    {
	return (v[0] <= (*this)[0] && v[1] <= (*this)[1] && v[2] <= (*this)[2]);
    }

    /**
     * componentwise greater than or equal to comparison
     */
    bool operator>=(vector3d<T> const& v) const
    {
	return (v[0] >= (*this)[0] && v[1] >= (*this)[1] && v[2] >= (*this)[2]);
    }

    /**
     * assignment by vector
     */
    vector3d<T>& operator=(vector3d<T> const& v)
    {
	(*this)[0] = v[0];
	(*this)[1] = v[1];
	(*this)[2] = v[2];
	return *this;
    }

    /**
     * assignment by scalar
     */
    vector3d<T>& operator=(T const& s)
    {
	(*this)[0] = s;
	(*this)[1] = s;
	(*this)[2] = s;
	return *this;
    }

    /**
     * assignment by componentwise vector addition
     */
    vector3d<T>& operator+=(vector3d<T> const& v)
    {
	(*this)[0] += v[0];
	(*this)[1] += v[1];
	(*this)[2] += v[2];
	return *this;
    }

    /**
     * assignment by componentwise vector subtraction
     */
    vector3d<T>& operator-=(vector3d<T> const& v)
    {
	(*this)[0] -= v[0];
	(*this)[1] -= v[1];
	(*this)[2] -= v[2];
	return *this;
    }

    /**
     * assignment by scalar multiplication
     */
    vector3d<T>& operator*=(T const& s)
    {
	(*this)[0] *= s;
	(*this)[1] *= s;
	(*this)[2] *= s;
	return *this;
    }

    /**
     * assignment by scalar division
     */
    vector3d<T>& operator/=(T const& s)
    {
	(*this)[0] /= s;
	(*this)[1] /= s;
	(*this)[2] /= s;
	return *this;
    }

    /**
     * componentwise vector addition
     */
    friend vector3d<T> operator+(vector3d<T> v, vector3d<T> const& w)
    {
	v[0] += w[0];
	v[1] += w[1];
	v[2] += w[2];
	return v;
    }

    /**
     * componentwise vector subtraction
     */
    friend vector3d<T> operator-(vector3d<T> v, vector3d<T> const& w)
    {
	v[0] -= w[0];
	v[1] -= w[1];
	v[2] -= w[2];
	return v;
    }

    /**
     * scalar product
     */
    T operator*(vector3d<T> const& v) const
    {
	return (*this)[0] * v[0] + (*this)[1] * v[1] + (*this)[2] * v[2];
    }

    /**
     * scalar multiplication
     */
    friend vector3d<T> operator*(vector3d<T> v, T const& s)
    {
	v[0] *= s;
	v[1] *= s;
	v[2] *= s;
	return v;
    }

    /**
     * scalar multiplication
     */
    friend vector3d<T> operator*(T const& s, vector3d<T> v)
    {
	v[0] *= s;
	v[1] *= s;
	v[2] *= s;
	return v;
    }

    /**
     * scalar division
     */
    friend vector3d<T> operator/(vector3d<T> v, T const& s)
    {
	v[0] /= s;
	v[1] /= s;
	v[2] /= s;
	return v;
    }

    /**
     * write vector components to output stream
     */
    friend std::ostream& operator<<(std::ostream& os, vector3d<T> const& v)
    {
	os << v[0] << "\t" << v[1] << "\t" << v[2];
	return os;
    }

    /**
     * read vector components from input stream
     */
    friend std::istream& operator>>(std::istream& is, vector3d<T>& v)
    {
	is >> v[0] >> v[1] >> v[2];
	return is;
    }

    /**
     * translate from vector to XDR vector
     */
    friend xdr::ostream& operator<<(xdr::ostream& xdrs, vector3d<T> const& v)
    {
	xdrs << v[0] << v[1] << v[2];
	return xdrs;
    }

    /**
     * translate to vector from XDR vector
     */
    friend xdr::istream& operator>>(xdr::istream& xdrs, vector3d<T>& v)
    {
	xdrs >> v[0] >> v[1] >> v[2];
	return xdrs;
    }
};

/**
 * componentwise round to nearest integer
 */
template <typename T>
vector3d<T> rint(vector3d<T> v);

template <>
vector3d<float> rint(vector3d<float> v)
{
    v[0] = rintf(v[0]);
    v[1] = rintf(v[1]);
    v[2] = rintf(v[2]);
    return v;
}

template <>
vector3d<double> rint(vector3d<double> v)
{
    v[0] = rint(v[0]);
    v[1] = rint(v[1]);
    v[2] = rint(v[2]);
    return v;
}

/**
 * componentwise round to nearest integer, away from zero
 */
template <typename T>
vector3d<T> round(vector3d<T> v);

template <>
vector3d<float> round(vector3d<float> v)
{
    v[0] = roundf(v[0]);
    v[1] = roundf(v[1]);
    v[2] = roundf(v[2]);
    return v;
}

template <>
vector3d<double> round(vector3d<double> v)
{
    v[0] = round(v[0]);
    v[1] = round(v[1]);
    v[2] = round(v[2]);
    return v;
}

/**
 * componentwise round to nearest integer not greater than argument
 */
template <typename T>
vector3d<T> floor(vector3d<T> v)
{
    v[0] = std::floor(v[0]);
    v[1] = std::floor(v[1]);
    v[2] = std::floor(v[2]);
    return v;
}

/**
 * componentwise round to nearest integer not less argument
 */
template <typename T>
vector3d<T> ceil(vector3d<T> v)
{
    v[0] = std::ceil(v[0]);
    v[1] = std::ceil(v[1]);
    v[2] = std::ceil(v[2]);
    return v;
}

/**
 * componentwise square root function
 */
template <typename T>
vector3d<T> sqrt(vector3d<T> v)
{
    v[0] = std::sqrt(v[0]);
    v[1] = std::sqrt(v[1]);
    v[2] = std::sqrt(v[2]);
    return v;
}

/**
 * componentwise cos function
 */
template <typename T>
vector3d<T> cos(vector3d<T> v)
{
    v[0] = std::cos(v[0]);
    v[1] = std::cos(v[1]);
    v[2] = std::cos(v[2]);
    return v;
}

/**
 * componentwise sin function
 */
template <typename T>
vector3d<T> sin(vector3d<T> v)
{
    v[0] = std::sin(v[0]);
    v[1] = std::sin(v[1]);
    v[2] = std::sin(v[2]);
    return v;
}

#endif /* ! MDSIM_VECTOR3D_HPP */
