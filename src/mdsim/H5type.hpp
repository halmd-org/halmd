/* HDF5 native data type translation
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

#include <H5Cpp.h>

template <typename T>
struct H5type;

template <>
struct H5type<int8_t> : public H5::DataType
{
    H5type() : H5::DataType(H5::PredType::NATIVE_INT8) {}
};

template <>
struct H5type<uint8_t> : public H5::DataType
{
    H5type() : H5::DataType(H5::PredType::NATIVE_UINT8) {}
};

template <>
struct H5type<int16_t> : public H5::DataType
{
    H5type() : H5::DataType(H5::PredType::NATIVE_INT16) {}
};

template <>
struct H5type<uint16_t> : public H5::DataType
{
    H5type() : H5::DataType(H5::PredType::NATIVE_UINT16) {}
};

template <>
struct H5type<int32_t> : public H5::DataType
{
    H5type() : H5::DataType(H5::PredType::NATIVE_INT32) {}
};

template <>
struct H5type<uint32_t> : public H5::DataType
{
    H5type() : H5::DataType(H5::PredType::NATIVE_UINT32) {}
};

template <>
struct H5type<int64_t> : public H5::DataType
{
    H5type() : H5::DataType(H5::PredType::NATIVE_INT64) {}
};

template <>
struct H5type<uint64_t> : public H5::DataType
{
    H5type() : H5::DataType(H5::PredType::NATIVE_UINT64) {}
};

template <>
struct H5type<float> : public H5::DataType
{
    H5type() : H5::DataType(H5::PredType::NATIVE_FLOAT) {}
};

template <>
struct H5type<double> : public H5::DataType
{
    H5type() : H5::DataType(H5::PredType::NATIVE_DOUBLE) {}
};
