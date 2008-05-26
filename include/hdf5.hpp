/* HDF5 C++ API extensions
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

#ifndef H5EXT_HPP
#define H5EXT_HPP

#include <H5Cpp.h>
#include <stdint.h>
#include <string.h>
#include <string>


/**
 * HDF5 C++ API extensions
 */
namespace H5ext
{

/**
 * HDF5 native datatype translation
 */
template <typename T>
struct NativeType;

template <>
struct NativeType<int8_t> : public H5::DataType
{
    NativeType() : H5::DataType(H5::PredType::NATIVE_INT8) {}
};

template <>
struct NativeType<uint8_t> : public H5::DataType
{
    NativeType() : H5::DataType(H5::PredType::NATIVE_UINT8) {}
};

template <>
struct NativeType<int16_t> : public H5::DataType
{
    NativeType() : H5::DataType(H5::PredType::NATIVE_INT16) {}
};

template <>
struct NativeType<uint16_t> : public H5::DataType
{
    NativeType() : H5::DataType(H5::PredType::NATIVE_UINT16) {}
};

template <>
struct NativeType<int32_t> : public H5::DataType
{
    NativeType() : H5::DataType(H5::PredType::NATIVE_INT32) {}
};

template <>
struct NativeType<uint32_t> : public H5::DataType
{
    NativeType() : H5::DataType(H5::PredType::NATIVE_UINT32) {}
};

template <>
struct NativeType<int64_t> : public H5::DataType
{
    NativeType() : H5::DataType(H5::PredType::NATIVE_INT64) {}
};

template <>
struct NativeType<uint64_t> : public H5::DataType
{
    NativeType() : H5::DataType(H5::PredType::NATIVE_UINT64) {}
};

template <>
struct NativeType<float> : public H5::DataType
{
    NativeType() : H5::DataType(H5::PredType::NATIVE_FLOAT) {}
};

template <>
struct NativeType<double> : public H5::DataType
{
    NativeType() : H5::DataType(H5::PredType::NATIVE_DOUBLE) {}
};


/**
 * HDF5 attribute
 */
class Attribute : public H5::Attribute
{
public:
    Attribute(H5::Group const& node, std::string const& name) : node_(node), name_(name) {}
    template <typename T> Attribute& operator=(T const& value);
    template <typename T> Attribute& operator=(T const* value);
    template <typename T> T as();

private:
    /** group which attribute belongs to */
    H5::Group node_;
    /** attribute name */
    std::string name_;
};

/*
 * create and write native attribute
 */
template <typename T>
Attribute& Attribute::operator=(T const& value)
{
    NativeType<T> tid;
    H5::Attribute::operator=(node_.createAttribute(name_, tid, H5S_SCALAR));
    H5::Attribute::write(tid, &value);
    return *this;
}

/**
 * create and write string attribute
 */
template <>
Attribute& Attribute::operator=(std::string const& value)
{
    H5::StrType tid(H5::PredType::C_S1, value.length());
    H5::Attribute::operator=(node_.createAttribute(name_, tid, H5S_SCALAR));
    H5::Attribute::write(tid, value.c_str());
    return *this;
}

/**
 * create and write string attribute
 */
template <>
Attribute& Attribute::operator=(char const* value)
{
    H5::StrType tid(H5::PredType::C_S1, strlen(value));
    H5::Attribute::operator=(node_.createAttribute(name_, tid, H5S_SCALAR));
    H5::Attribute::write(tid, value);
    return *this;
}

/**
 * read native attribute
 */
template <typename T>
T Attribute::as()
{
    NativeType<T> tid;
    T value;
    H5::Attribute::operator=(node_.openAttribute(name_));
    H5::Attribute::read(tid, &value);
    return value;
}

/**
 * read string attribute
 */
template <>
std::string Attribute::as()
{
    H5::Attribute::operator=(node_.openAttribute(name_));

    //
    // Instead of using the returned datatype of the attribute,
    // we explicitly define a string datatype of approriate size
    // to let the HDF5 library handle a possible type mismatch.
    // 
    H5::StrType tid(H5::PredType::C_S1, H5::Attribute::getDataType().getSize());

    std::string value(tid.getSize(), '\0');
    H5::Attribute::read(tid, &(*value.begin()));
    return value;
}


/**
 * HDF5 group
 */
class Group : public H5::Group
{
public:
    Group() {}
    Group(H5::Group const& node) : H5::Group(node) {}
    Attribute operator[](char const* name);
};

/**
 * create and write an attribute
 */
Attribute Group::operator[](char const* name)
{
    return Attribute(*this, name);
}

} // namespace H5ext

#endif /* ! H5EXT_HPP */
