/* HDF5 C++ extensions
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

#ifndef LJGPU_UTIL_H5XX_HPP
#define LJGPU_UTIL_H5XX_HPP

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
    attribute(H5::Group const& node, std::string const& name)
	: m_node(node), m_name(name) {}
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

template <typename Exception>
class no_autoprint : public Exception
{
public:
    no_autoprint()
    {
	Exception::getAutoPrint(func, &client_data);
	Exception::dontPrint();
    }

    ~no_autoprint()
    {
	Exception::setAutoPrint(func, client_data);
    }

private:
    H5E_auto_t func;
    void* client_data;
};

#define H5XX_NO_AUTO_PRINT(exception) H5xx::no_autoprint<exception> __no_autoprint;

} // namespace H5xx

#endif /* ! LJGPU_UTIL_H5XX_HPP */
