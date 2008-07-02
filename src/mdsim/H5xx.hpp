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

#ifndef MDSIM_H5XX_HPP
#define MDSIM_H5XX_HPP

#include <H5Cpp.h>

namespace H5xx
{

/*
 * C type to HDF5 native data type translation
 */
template <typename T>
struct ctype
{
    static H5::PredType const& type;
};

/**
 * HDF5 attribute
 */
class attribute
{
public:
    attribute(H5::Group const& node, std::string const& name) : m_node(node), m_name(name)
    {
	// turns off the automatic error printing from the HDF5 library
	H5::AttributeIException::dontPrint();
    }

    template <typename T> attribute& operator=(T const& value);
    template <typename T> attribute& operator=(T const* value);
    template <typename T> T as();

private:
    /** group which attribute belongs to */
    H5::Group m_node;
    /** attribute name */
    std::string m_name;
};

/**
 * HDF5 group
 */
class group : public H5::Group
{
public:
    group() {}
    group(H5::Group const& node) : H5::Group(node) {}

    /**
     * returns existing or creates attribute
     */
    attribute operator[](char const* name) const
    {
	return attribute(*this, name);
    }
};

} // namespace H5xx

#endif /* ! H5XX_HPP */
