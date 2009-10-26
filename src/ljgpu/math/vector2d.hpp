/* 2-dimensional floating-point vector
 *
 * Copyright © 2008-2009  Peter Colberg
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

#ifndef LJGPU_MATH_VECTOR2D_HPP
#define LJGPU_MATH_VECTOR2D_HPP

#include <boost/array.hpp>
#include <cmath>
#include <iostream>
#ifdef WITH_CUDA
# include <cuda_runtime.h>
#endif

// overloaded math functions must be defined in global namespace

template <typename T, unsigned int dimension>
class vector;

/**
 * 2-dimensional floating-point vector
 */
template <typename T>
class vector<T, 2> : public boost::array<T, 2>
{
public:
    typedef T value_type;

public:
    vector() {}

    /**
     * initialization by vector
     */
    template <typename U>
    vector(vector<U, 2> const& v)
    {
        (*this)[0] = v[0];
        (*this)[1] = v[1];
    }

    /**
     * initialization by GPU floating-point vector
     */
#ifdef WITH_CUDA
    vector(float2 const& v)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
    }

    vector(float3 const& v)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
    }

    vector(float4 const& v)
    {
        (*this)[0] = v.x;
        (*this)[1] = v.y;
    }
#endif /* WITH_CUDA */

    /**
     * initialization by scalar
     */
    vector(T const& s)
    {
        (*this)[0] = s;
        (*this)[1] = s;
    }

    /**
     * initialization by scalar components
     */
    vector(T const& x, T const& y)
    {
        (*this)[0] = x;
        (*this)[1] = y;
    }

    /**
     * convert to GPU floating-point type
     */
#ifdef WITH_CUDA
    operator float2() const
    {
        return make_float2((*this)[0], (*this)[1]);
    }

    operator float3() const
    {
        return make_float3((*this)[0], (*this)[1], 0);
    }

    operator float4() const
    {
        return make_float4((*this)[0], (*this)[1], 0, 0);
    }
#endif /* WITH_CUDA */

    /**
     * equality comparison
     */
    bool operator==(vector<T, 2> const& v) const
    {
        return (v[0] == (*this)[0] && v[1] == (*this)[1]);
    }

    /**
     * inequality comparison
     */
    bool operator!=(vector<T, 2> const& v) const
    {
        return (v[0] != (*this)[0] || v[1] != (*this)[1]);
    }

    /**
     * componentwise less than comparison
     */
    bool operator<(vector<T, 2> const& v) const
    {
        return (v[0] < (*this)[0] && v[1] < (*this)[1]);
    }

    /**
     * componentwise greater than comparison
     */
    bool operator>(vector<T, 2> const& v) const
    {
        return (v[0] > (*this)[0] && v[1] > (*this)[1]);
    }

    /**
     * componentwise less than or equal to comparison
     */
    bool operator<=(vector<T, 2> const& v) const
    {
        return (v[0] <= (*this)[0] && v[1] <= (*this)[1]);
    }

    /**
     * componentwise greater than or equal to comparison
     */
    bool operator>=(vector<T, 2> const& v) const
    {
        return (v[0] >= (*this)[0] && v[1] >= (*this)[1]);
    }

    /**
     * assignment by vector
     */
    vector<T, 2>& operator=(vector<T, 2> const& v)
    {
        (*this)[0] = v[0];
        (*this)[1] = v[1];
        return *this;
    }

    /**
     * assignment by scalar
     */
    vector<T, 2>& operator=(T const& s)
    {
        (*this)[0] = s;
        (*this)[1] = s;
        return *this;
    }

    /**
     * assignment by componentwise vector addition
     */
    vector<T, 2>& operator+=(vector<T, 2> const& v)
    {
        (*this)[0] += v[0];
        (*this)[1] += v[1];
        return *this;
    }

    /**
     * assignment by componentwise vector subtraction
     */
    vector<T, 2>& operator-=(vector<T, 2> const& v)
    {
        (*this)[0] -= v[0];
        (*this)[1] -= v[1];
        return *this;
    }

    /**
     * assignment by scalar multiplication
     */
    vector<T, 2>& operator*=(T const& s)
    {
        (*this)[0] *= s;
        (*this)[1] *= s;
        return *this;
    }

    /**
     * assignment by scalar division
     */
    vector<T, 2>& operator/=(T const& s)
    {
        (*this)[0] /= s;
        (*this)[1] /= s;
        return *this;
    }

    /**
     * componentwise vector addition
     */
    friend vector<T, 2> operator+(vector<T, 2> v, vector<T, 2> const& w)
    {
        v[0] += w[0];
        v[1] += w[1];
        return v;
    }

    /**
     * componentwise vector subtraction
     */
    friend vector<T, 2> operator-(vector<T, 2> v, vector<T, 2> const& w)
    {
        v[0] -= w[0];
        v[1] -= w[1];
        return v;
    }

    /**
     * scalar product
     */
    T operator*(vector<T, 2> const& v) const
    {
        return (*this)[0] * v[0] + (*this)[1] * v[1];
    }

    /**
     * scalar multiplication
     */
    friend vector<T, 2> operator*(vector<T, 2> v, T const& s)
    {
        v[0] *= s;
        v[1] *= s;
        return v;
    }

    /**
     * scalar multiplication
     */
    friend vector<T, 2> operator*(T const& s, vector<T, 2> v)
    {
        v[0] *= s;
        v[1] *= s;
        return v;
    }

    /**
     * scalar division
     */
    friend vector<T, 2> operator/(vector<T, 2> v, T const& s)
    {
        v[0] /= s;
        v[1] /= s;
        return v;
    }

    /**
     * write vector components to output stream
     */
    friend std::ostream& operator<<(std::ostream& os, vector<T, 2> const& v)
    {
        os << v[0] << "\t" << v[1];
        return os;
    }

    /**
     * read vector components from input stream
     */
    friend std::istream& operator>>(std::istream& is, vector<T, 2>& v)
    {
        is >> v[0] >> v[1];
        return is;
    }
};

/**
 * componentwise round to nearest integer
 */
template <typename T>
vector<T, 2> rint(vector<T, 2> v);

template <>
inline vector<float, 2> rint(vector<float, 2> v)
{
    v[0] = ::rintf(v[0]);
    v[1] = ::rintf(v[1]);
    return v;
}

template <>
inline vector<double, 2> rint(vector<double, 2> v)
{
    v[0] = ::rint(v[0]);
    v[1] = ::rint(v[1]);
    return v;
}

/**
 * componentwise round to nearest integer, away from zero
 */
template <typename T>
vector<T, 2> round(vector<T, 2> v);

template <>
inline vector<float, 2> round(vector<float, 2> v)
{
    v[0] = ::roundf(v[0]);
    v[1] = ::roundf(v[1]);
    return v;
}

template <>
inline vector<double, 2> round(vector<double, 2> v)
{
    v[0] = ::round(v[0]);
    v[1] = ::round(v[1]);
    return v;
}

/**
 * componentwise round to nearest integer not greater than argument
 */
template <typename T>
inline vector<T, 2> floor(vector<T, 2> v)
{
    v[0] = std::floor(v[0]);
    v[1] = std::floor(v[1]);
    return v;
}

/**
 * componentwise round to nearest integer not less argument
 */
template <typename T>
inline vector<T, 2> ceil(vector<T, 2> v)
{
    v[0] = std::ceil(v[0]);
    v[1] = std::ceil(v[1]);
    return v;
}

/**
 * componentwise square root function
 */
template <typename T>
inline vector<T, 2> sqrt(vector<T, 2> v)
{
    v[0] = std::sqrt(v[0]);
    v[1] = std::sqrt(v[1]);
    return v;
}

/**
 * componentwise cos function
 */
template <typename T>
inline vector<T, 2> cos(vector<T, 2> v)
{
    v[0] = std::cos(v[0]);
    v[1] = std::cos(v[1]);
    return v;
}

/**
 * componentwise sin function
 */
template <typename T>
inline vector<T, 2> sin(vector<T, 2> v)
{
    v[0] = std::sin(v[0]);
    v[1] = std::sin(v[1]);
    return v;
}

#endif /* ! LJGPU_MATH_VECTOR2D_HPP */
