/* HDF5 parameter group
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

#ifndef MDSIM_H5PARAM_HPP
#define MDSIM_H5PARAM_HPP

#include <H5Cpp.h>


namespace mdsim
{

/**
 * HDF5 parameter group
 */
class H5param
{
public:
    /** initialize HDF5 parameter group */
    H5param(H5::Group param);
    /** write parameters of visitor to HDF5 parameter group */
    template <typename T>
    H5param& operator<<(T const& visitor);

private:
    H5::Group param;
};

/**
 * write parameters of visitor to HDF5 parameter group
 */
template <typename T>
H5param& H5param::operator<<(T const& visitor)
{
    visitor.attrs(param);
    return *this;
}

} // namespace mdsim

#endif /* ! MDSIM_H5PARAM_HPP */
