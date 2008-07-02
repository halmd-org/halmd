/* HDF5 C++ extensions
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

#include <stdint.h>
#include "H5xx.hpp"

namespace H5xx
{

/*
 * create and write C type attribute
 */
template <typename T>
attribute& attribute::operator=(T const& value)
{
    H5::Attribute attr;
    try {
	attr = m_node.openAttribute(m_name);
    }
    catch (H5::AttributeIException const&) {
	attr = m_node.createAttribute(m_name, ctype<T>::type, H5S_SCALAR);
    }
    attr.write(ctype<T>::type, &value);
    return *this;
}

/**
 * create and write string attribute
 */
template <>
attribute& attribute::operator=(std::string const& value)
{
    H5::StrType tid(H5::PredType::C_S1, 256);
    H5::Attribute attr;
    try {
	attr = m_node.openAttribute(m_name);
    }
    catch (H5::AttributeIException const&) {
	attr = m_node.createAttribute(m_name, tid, H5S_SCALAR);
    }
    attr.write(tid, value.c_str());
    return *this;
}

/**
 * create and write C string attribute
 */
template <>
attribute& attribute::operator=(char const* value)
{
    H5::StrType tid(H5::PredType::C_S1, 256);
    H5::Attribute attr;
    try {
	attr = m_node.openAttribute(m_name);
    }
    catch (H5::AttributeIException const&) {
	attr = m_node.createAttribute(m_name, tid, H5S_SCALAR);
    }
    attr.write(tid, value);
    return *this;
}

/**
 * read C type attribute
 */
template <typename T>
T attribute::as()
{
    H5::Attribute attr;
    try {
	attr = m_node.openAttribute(m_name);
    }
    catch (H5::AttributeIException const&) {
	throw;
    }
    T value;
    attr.read(ctype<T>::type, &value);
    return value;
}

/**
 * read string attribute
 */
template <>
std::string attribute::as()
{
    H5::Attribute attr;
    try {
	attr = m_node.openAttribute(m_name);
    }
    catch (H5::AttributeIException const&) {
	throw;
    }
    // fixed string length includes terminating NULL character
    char value[256];
    attr.read(H5::StrType(H5::PredType::C_S1, 256), value);
    return value;
}

template <> H5::PredType const& ctype<int8_t>::type = H5::PredType::NATIVE_INT8;
template <> H5::PredType const& ctype<uint8_t>::type = H5::PredType::NATIVE_UINT8;
template <> H5::PredType const& ctype<int16_t>::type = H5::PredType::NATIVE_INT16;
template <> H5::PredType const& ctype<uint16_t>::type = H5::PredType::NATIVE_UINT16;
template <> H5::PredType const& ctype<int32_t>::type = H5::PredType::NATIVE_INT32;
template <> H5::PredType const& ctype<uint32_t>::type = H5::PredType::NATIVE_UINT32;
template <> H5::PredType const& ctype<int64_t>::type = H5::PredType::NATIVE_INT64;
template <> H5::PredType const& ctype<uint64_t>::type = H5::PredType::NATIVE_UINT64;
template <> H5::PredType const& ctype<float>::type = H5::PredType::NATIVE_FLOAT;
template <> H5::PredType const& ctype<double>::type = H5::PredType::NATIVE_DOUBLE;

// explicit template instantiation
template attribute& attribute::operator=(int8_t const&);
template attribute& attribute::operator=(uint8_t const&);
template attribute& attribute::operator=(int16_t const&);
template attribute& attribute::operator=(uint16_t const&);
template attribute& attribute::operator=(int32_t const&);
template attribute& attribute::operator=(uint32_t const&);
template attribute& attribute::operator=(int64_t const&);
template attribute& attribute::operator=(uint64_t const&);
template attribute& attribute::operator=(float const&);
template attribute& attribute::operator=(double const&);
template int8_t attribute::as();
template uint8_t attribute::as();
template int16_t attribute::as();
template uint16_t attribute::as();
template int32_t attribute::as();
template uint32_t attribute::as();
template int64_t attribute::as();
template uint64_t attribute::as();
template float attribute::as();
template double attribute::as();

} // namespace H5xx
