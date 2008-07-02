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

#include <stdint.h>
#include "H5param.hpp"
#include "H5type.hpp"
#include "exception.hpp"
#include "version.h"


namespace mdsim
{

/**
 * initialize HDF5 parameter group
 */
H5param::H5param(H5::Group param) : param(param)
{
    H5::DataType tid(H5::StrType(H5::PredType::C_S1, 256));
    // write program info attributes
    H5::Group node(param.createGroup("program"));
    node.createAttribute("name", tid, H5S_SCALAR).write(tid, PROGRAM_NAME);
    node.createAttribute("version", tid, H5S_SCALAR).write(tid, PROGRAM_VERSION);
    node.createAttribute("variant", tid, H5S_SCALAR).write(tid, PROGRAM_VARIANT);
}

/**
 * create attribute in given HDF5 group
 */
template <typename T>
void H5param::attr(H5::Group const& node, char const* name, T value)
{
    try {
	H5type<T> tid;
	H5::Attribute attr(node.createAttribute(name, tid, H5S_SCALAR));
	attr.write(tid, &value);
    }
    catch (H5::Exception const& e) {
	throw exception("failed to write parameter to HDF5 output file");
    }
}

// explicit template instantiation
template void H5param::attr(H5::Group const&, char const*, int8_t);
template void H5param::attr(H5::Group const&, char const*, uint8_t);
template void H5param::attr(H5::Group const&, char const*, int16_t);
template void H5param::attr(H5::Group const&, char const*, uint16_t);
template void H5param::attr(H5::Group const&, char const*, int32_t);
template void H5param::attr(H5::Group const&, char const*, uint32_t);
template void H5param::attr(H5::Group const&, char const*, int64_t);
template void H5param::attr(H5::Group const&, char const*, uint64_t);
template void H5param::attr(H5::Group const&, char const*, float);
template void H5param::attr(H5::Group const&, char const*, double);

} // namespace mdsim
