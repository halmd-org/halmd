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
#include <cuda_wrapper.hpp>
#include <xdr/iostream.hpp>


/**
 * 3-dimensional floating-point vector
 */
template <typename T>
class vector3d
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
    vector3d(vector3d<T> const& v) : x(v.x), y(v.y), z(v.z)
    {
    }

    /**
     * initialization by vector
     */
    vector3d(float3 const& v) : x(v.x), y(v.y), z(v.z)
    {
    }

    /**
     * initialization by scalar
     */
    vector3d(T const& s) : x(s), y(s), z(s)
    {
    }

    /**
     * initialization by scalar components
     */
    vector3d(T const& x, T const& y, T const& z) : x(x), y(y), z(z)
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
     * equality comparison
     */
    bool operator==(vector3d<T> const& v) const
    {
        return (v.x == x && v.y == y && v.z == z);
    }
    
    /**
     * inequality comparison
     */
    bool operator!=(vector3d<T> const& v) const
    {
        return (v.x != x || v.y != y || v.z != z);
    }

    /**
     * componentwise less than comparison
     */
    bool operator<(vector3d<T> const& v) const
    {
	return (v.x < x && v.y < y && v.z < z);
    }

    /**
     * componentwise greater than comparison
     */
    bool operator>(vector3d<T> const& v) const
    {
	return (v.x > x && v.y > y && v.z > z);
    }

    /**
     * componentwise less than or equal to comparison
     */
    bool operator<=(vector3d<T> const& v) const
    {
	return (v.x <= x && v.y <= y && v.z <= z);
    }

    /**
     * componentwise greater than or equal to comparison
     */
    bool operator>=(vector3d<T> const& v) const
    {
	return (v.x >= x && v.y >= y && v.z >= z);
    }

    /**
     * assignment by vector
     */
    vector3d<T>& operator=(vector3d<T> const& v)
    {
	x = v.x;
	y = v.y;
	z = v.z;
	return *this;
    }

    /**
     * assignment by scalar
     */
    vector3d<T>& operator=(T const& s)
    {
	x = s;
	y = s;
	z = s;
	return *this;
    }

    /**
     * assignment by componentwise vector addition
     */
    vector3d<T>& operator+=(vector3d<T> const& v)
    {
	x += v.x;
	y += v.y;
	z += v.z;
	return *this;
    }

    /**
     * assignment by componentwise vector subtraction
     */
    vector3d<T>& operator-=(vector3d<T> const& v)
    {
	x -= v.x;
	y -= v.y;
	z -= v.z;
	return *this;
    }

    /**
     * assignment by scalar multiplication
     */
    vector3d<T>& operator*=(T const& s)
    {
	x *= s;
	y *= s;
	z *= s;
	return *this;
    }

    /**
     * assignment by scalar division
     */
    vector3d<T>& operator/=(T const& s)
    {
	x /= s;
	y /= s;
	z /= s;
	return *this;
    }

    /**
     * componentwise vector addition
     */
    friend vector3d<T> operator+(vector3d<T> v, vector3d<T> const& w)
    {
	v.x += w.x;
	v.y += w.y;
	v.z += w.z;
	return v;
    }

    /**
     * componentwise vector subtraction
     */
    friend vector3d<T> operator-(vector3d<T> v, vector3d<T> const& w)
    {
	v.x -= w.x;
	v.y -= w.y;
	v.z -= w.z;
	return v;
    }

    /**
     * scalar product
     */
    T operator*(vector3d<T> const& v) const
    {
	return x * v.x + y * v.y + z * v.z;
    }

    /**
     * scalar multiplication
     */
    friend vector3d<T> operator*(vector3d<T> v, T const& s)
    {
	v.x *= s;
	v.y *= s;
	v.z *= s;
	return v;
    }

    /**
     * scalar multiplication
     */
    friend vector3d<T> operator*(T const& s, vector3d<T> v)
    {
	v.x *= s;
	v.y *= s;
	v.z *= s;
	return v;
    }

    /**
     * scalar division
     */
    friend vector3d<T> operator/(vector3d<T> v, T const& s)
    {
	v.x /= s;
	v.y /= s;
	v.z /= s;
	return v;
    }

    /**
     * write vector components to output stream
     */
    friend std::ostream& operator<<(std::ostream& os, vector3d<T> const& v)
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

    /**
     * translate from vector to XDR vector
     */
    friend xdr::ostream& operator<<(xdr::ostream& xdrs, vector3d<T> const& v)
    {
	xdrs << v.x << v.y << v.z;
	return xdrs;
    }

    /**
     * translate to vector from XDR vector
     */
    friend xdr::istream& operator>>(xdr::istream& xdrs, vector3d<T>& v)
    {
	xdrs >> v.x >> v.y >> v.z;
	return xdrs;
    }

public:
    T x, y, z;
};


/**
 * returns device pointer to allocated device memory
 */
float3* cuda_cast(cuda::vector<vector3d<float> >& v)
{
    return reinterpret_cast<float3*>(v.data());
}

/**
 * returns device pointer to allocated device memory
 */
float3 const* cuda_cast(cuda::vector<vector3d<float> > const& v)
{
    return reinterpret_cast<float3 const*>(v.data());
}


/**
 * componentwise round to nearest integer
 */
template <typename T>
vector3d<T> rint(vector3d<T> v);

template <>
vector3d<float> rint(vector3d<float> v)
{
    v.x = rintf(v.x);
    v.y = rintf(v.y);
    v.z = rintf(v.z);
    return v;
}

template <>
vector3d<double> rint(vector3d<double> v)
{
    v.x = rint(v.x);
    v.y = rint(v.y);
    v.z = rint(v.z);
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
    v.x = roundf(v.x);
    v.y = roundf(v.y);
    v.z = roundf(v.z);
    return v;
}

template <>
vector3d<double> round(vector3d<double> v)
{
    v.x = round(v.x);
    v.y = round(v.y);
    v.z = round(v.z);
    return v;
}


/**
 * componentwise round to nearest integer not greater than argument
 */
template <typename T>
vector3d<T> floor(vector3d<T> v);

template <>
vector3d<float> floor(vector3d<float> v)
{   
    v.x = floorf(v.x);
    v.y = floorf(v.y);
    v.z = floorf(v.z);
    return v;
}

template <>
vector3d<double> floor(vector3d<double> v)
{   
    v.x = floor(v.x);
    v.y = floor(v.y);
    v.z = floor(v.z);
    return v;
}


/**
 * componentwise round to nearest integer not less argument
 */
template <typename T>
vector3d<T> ceil(vector3d<T> v);

template <>
vector3d<float> ceil(vector3d<float> v)
{   
    v.x = ceilf(v.x);
    v.y = ceilf(v.y);
    v.z = ceilf(v.z);
    return v;
}

template <>
vector3d<double> ceil(vector3d<double> v)
{   
    v.x = ceil(v.x);
    v.y = ceil(v.y);
    v.z = ceil(v.z);
    return v;
}


/**
 * componentwise square root function
 */
template <typename T>
vector3d<T> sqrt(vector3d<T> v);

template <>
vector3d<float> sqrt(vector3d<float> v)
{   
    v.x = sqrtf(v.x);
    v.y = sqrtf(v.y);
    v.z = sqrtf(v.z);
    return v;
}

template <>
vector3d<double> sqrt(vector3d<double> v)
{   
    v.x = sqrt(v.x);
    v.y = sqrt(v.y);
    v.z = sqrt(v.z);
    return v;
}

#endif /* ! MDSIM_VECTOR3D_HPP */
