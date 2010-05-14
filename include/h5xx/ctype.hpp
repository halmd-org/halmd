/*
 * Copyright © 2010  Peter Colberg
 *
 * This file is part of HALMD.
 *
 * HALMD is free software: you can redistribute it and/or modify
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

#ifndef H5XX_CTYPE_HPP
#define H5XX_CTYPE_HPP

#include <h5xx/compat.hpp>

namespace h5xx
{

/**
 * HDF5 attribute
 */
template <typename T>
struct ctype;

template <>
struct ctype<float>
{
    operator hid_t() const
    {
        return H5T_NATIVE_FLOAT;
    }
};

template <>
struct ctype<double>
{
    operator hid_t() const
    {
        return H5T_NATIVE_DOUBLE;
    }
};

template <>
struct ctype<char>
{
    operator hid_t() const
    {
        return H5T_NATIVE_CHAR;
    }
};

template <>
struct ctype<signed char>
{
    operator hid_t() const
    {
        return H5T_NATIVE_SCHAR;
    }
};

template <>
struct ctype<unsigned char>
{
    operator hid_t() const
    {
        return H5T_NATIVE_UCHAR;
    }
};

template <>
struct ctype<short>
{
    operator hid_t() const
    {
        return H5T_NATIVE_SHORT;
    }
};

template <>
struct ctype<unsigned short>
{
    operator hid_t() const
    {
        return H5T_NATIVE_USHORT;
    }
};

template <>
struct ctype<int>
{
    operator hid_t() const
    {
        return H5T_NATIVE_INT;
    }
};

template <>
struct ctype<unsigned int>
{
    operator hid_t() const
    {
        return H5T_NATIVE_UINT;
    }
};

template <>
struct ctype<long>
{
    operator hid_t() const
    {
        return H5T_NATIVE_LONG;
    }
};

template <>
struct ctype<unsigned long>
{
    operator hid_t() const
    {
        return H5T_NATIVE_ULONG;
    }
};

template <>
struct ctype<long long>
{
    operator hid_t() const
    {
        return H5T_NATIVE_LLONG;
    }
};

template <>
struct ctype<unsigned long long>
{
    operator hid_t() const
    {
        return H5T_NATIVE_ULLONG;
    }
};

} // namespace h5xx

#endif /* ! H5XX_CTYPE_HPP */
